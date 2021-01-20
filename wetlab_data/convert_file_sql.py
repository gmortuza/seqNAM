import math
import sqlite3 as sl

from Bio import SeqIO
import multiprocessing
from functools import partial
from collections import defaultdict
import pickle

connection = sl.connect("sequences.db")
cursor = connection.cursor()


def round_nearest_10(number):
    if number % 10 < 5:
        # round down
        return int(number/10) * 10
    else:
        # round up
        return int((number + 10) / 10) * 10


def analyze_seq_length():
    # Right now this function uses read which have the length of the sequences.
    # TODO: Combine this method with rest of the database query

    counter_map = {
        "240": "235 - 244",
        "230": "225 - 234",
        "220": "215 - 224",
        "210": "205 - 214",
        "260": "255 - 264",
        "270": "265 - 274",
        "280": "275 - 284",
        "290": "285 - 294",
        "less": "0 - 204",
        "greater": "295 - "
    }
    counter = defaultdict(lambda: 0)
    total = 0
    with open("test.txt", "r") as file:
        for line in file:
            length_rounded = int(math.ceil(int(line)/10.0)) * 10
            if length_rounded == 250:
                counter[line.rstrip()] += 1
            elif length_rounded < 210:
                counter[counter_map["less"]] += 1
            elif length_rounded > 290:
                counter[counter_map["greater"]] += 1
            else:
                counter[counter_map[str(length_rounded)]] += 1
            total += 1

    # sort based on value
    counter = dict(sorted(counter.items(), key=lambda item: item[1]))

    length_counter = dict()
    length_counter["raw"] = counter
    length_counter["percentage"] = {k: v * 100 / total for k, v in counter.items()}

    with open("length_counter.pkl", "wb") as counter_file:
        pickle.dump(length_counter, counter_file)


def get_core_seq(seq, forward_primer="ACATCCAACACTCTACGCCC", backward_primer="GAATAGGAGCCGCAACACAC",
                 based_on_primer_length=False):
    if based_on_primer_length:
        return seq[len(forward_primer): - len(backward_primer)]
    else:
        # we will remove the sequences based on the content of the primer
        # the last 5 parts of our primer is CGCCC so if we find this element in the sequences that means that's the end
        # of the primer
        try:
            first_primer_end_index = seq.index("CCC") + 3
            seq = seq[first_primer_end_index:]
            # remove the last primer
            backward_primer_begin_index = seq.rindex("GAAT")
            seq = seq[:backward_primer_begin_index]
        except ValueError:
            return -1
        return seq


def create_table():
    # Create table if not exists
    # length --> length of the sequence with primer
    # r1_count --> # times this sequence appeared in read 1
    # r2_count --> # times this sequence appeared in read 2
    # gt_including_primer --> Id of the ground truth sequence. if this value is -1 then it means sequence contains err
    # gt_pc --> remove primer based on content and see if the sequences is correct
    # gt_pl --> remove primer based on length and see if the sequence is correct
    # probable_gt --> Ground truth after fixing the error. It is determined by measuring the minimum edit distance
    # probable_err --> probable error in the whole sequences
    # probable_err_in_primer --> probable error in the primer region
    # probable_err_pos --> position of the probable error
    # sequence --> The sequence that we are dealing with
    command = """
    CREATE TABLE IF NOT EXISTS sequences
    (
    	id INTEGER NOT NULL
    		CONSTRAINT sequences_pk
    			PRIMARY KEY AUTOINCREMENT,
    	length INTEGER NOT NULL,
    	r1_count INTEGER DEFAULT 0, 
    	r2_count INTEGER DEFAULT 0,
    	gt_including_primer INTEGER DEFAULT -1,
    	gt_pc INTEGER DEFAULT -1,
    	gt_pl INTEGER DEFAULT -1,
    	probable_gt INTEGER DEFAULT 0,
    	probable_err INTEGER DEFAULT -1,
    	probable_err_in_primer INTEGER DEFAULT -1,
    	probable_err_pos TEXT DEFAULT '',
    	sequence TEXT NOT NULL UNIQUE 
    );
    """
    connection.execute(command)
    # create that unique index
    connection.execute("CREATE UNIQUE INDEX IF NOT EXISTS sequences_id_uindex on sequences (id);")

    # create table for ground truth sequences
    command = """
    CREATE TABLE IF NOT EXISTS gt_sequences
    (
        id INTEGER NOT NULL CONSTRAINT gt_sequences_pk PRIMARY KEY AUTOINCREMENT,
        sequence_with_primer TEXT NOT NULL UNIQUE,
        sequence_without_primer TEXT NOT NULL UNIQUE
    );
    """
    connection.execute(command)
    # create unique index
    connection.execute("CREATE UNIQUE INDEX IF NOT EXISTS gt_sequences_id_uindex on gt_sequences (id);")


# insert the data fro ground truth
def insert_ground_truth(file):
    with open(file, "r") as csv_file:
        # first line is header
        next(csv_file)
        for line in csv_file:
            seq_with_primer = line.split(",")[-1].rstrip()
            seq_without_primer = seq_with_primer[20:-20]
            # insert these line into the database
            command = f"INSERT INTO gt_sequences (sequence_with_primer, sequence_without_primer) VALUES " \
                      f"('{seq_with_primer}', '{seq_without_primer}')"
            connection.execute(command)
        # save the database in disk
        connection.commit()


# Create table
create_table()
# Populate ground truth sequences
insert_ground_truth("gt_sequences.csv")


# store ground truth information into hashtable to reduce database query
def get_ground_truth():
    cursor.execute("SELECT id, sequence_without_primer, sequence_with_primer FROM gt_sequences")
    results = cursor.fetchall()
    without_primer = {}
    with_primer = {}
    for result in results:
        without_primer[result[1]] = result[0]  # seq: id
        with_primer[result[2]] = result[0]  # seq: id
    return with_primer, without_primer


gt_seq_with_primer, gt_seq_without_primer = get_ground_truth()


# insert fastq file content
def insert_wetlab_data(file, read=1):
    inserted_so_far = 1
    for seq_record in SeqIO.parse(file, "fastq"):
        seq = str(seq_record.seq)
        # if sequences already exists then update it's count only.
        # as we have already analyzed this sequence
        cursor.execute(f"SELECT * FROM sequences WHERE sequence='{seq}'")
        existing_seq = cursor.fetchall()
        if existing_seq:
            # We have already processed this sequence so just increment this counter
            if read == 1:
                # increment count for read 1
                cursor.execute(
                    f"UPDATE sequences SET r1_count = {existing_seq[0][2] + 1} WHERE id={existing_seq[0][0]}")
            elif read == 2:
                # increment count for read 2
                cursor.execute(
                    f"UPDATE sequences SET r2_count = {existing_seq[0][3] + 1} WHERE id={existing_seq[0][0]}")
        else:
            values = {
                'length': len(seq),
                'r1_count': 1 if read == 1 else 0,
                'r2_count': 1 if read == 2 else 0,
                'sequence': seq
            }
            # This is a new sequence that we need to analyze
            # we will search this sequences in gt_sequences
            gt_id = gt_seq_with_primer.get(seq, -1)
            if gt_id:
                # Insert this sequence this doesn't have any error
                values['gt_including_primer'] = gt_id
                values['gt_pc'] = gt_id
                values['gt_pl'] = gt_id
                values['probable_gt'] = gt_id
                values['probable_err'] = 0
                values['probable_err_in_primer'] = 0
                values['probable_err_pos'] = ''
                cmd = f"INSERT INTO sequences (length, r1_count, r2_count, gt_including_primer, gt_pc, gt_pl, probable_gt, " \
                      f"probable_err, probable_err_in_primer, probable_err_pos, sequence) VALUES (:length, :r1_count, :r2_count, " \
                      f":gt_including_primer, :gt_pc,:gt_pl, :probable_gt, :probable_err, :probable_err_in_primer, " \
                      f":probable_err_pos, :sequence)"
            else:
                values['gt_including_primer'] = -1  # sequence didn't match

                # check by removing the primer based on length
                # check if this sequence exists
                gt_id_length = gt_seq_without_primer.get(get_core_seq(seq, based_on_primer_length=True), -1)
                if gt_id_length:
                    # sequence exists
                    values['gt_pl'] = gt_id_length
                    values['probable_gt'] = gt_id_length
                else:
                    values['gt_pl'] = -1

                gt_id_content = gt_seq_without_primer.get(get_core_seq(seq, based_on_primer_length=False), -1)
                if gt_id_content:
                    # sequence exists
                    values['gt_pc'] = gt_id_content
                    values['probable_gt'] = gt_id_content
                else:
                    values['gt_pc'] = -1

                cmd = f"INSERT INTO sequences (length, r1_count, r2_count, gt_including_primer, gt_pc, gt_pl, " \
                      f"sequence) VALUES (:length, :r1_count, :r2_count, :gt_including_primer, :gt_pc,:gt_pl, :sequence)"
            cursor.execute(cmd, values)
        # Don't commit to the database offen
        # Commit to the database after handling 10,000 sequences
        if inserted_so_far % 10000 == 0:
            print(f"Read: {read}: sequence inserted: {inserted_so_far}")
            connection.commit()


def analyze_wetlab_data():
    # the first parameter of get_edit_distance method is (id, seq)
    # we don't require any id in the primer so we are passing -1
    # TODO: add this in analyzing sequence part
    # _, values['probable_err_in_primer'] = get_edit_distance((-1, "ACATCCAACACTCTACGCCCGAATAGGAGCCGCAACACAC"), seq[:20] + seq[-20:])
    # if we didn't find the sequence now need to perform exhaustive search on all the sequences
    # to get the minimum edit distance sequences
    # if not gt_id_content and not gt_id_length:
    #     # Minimum edit distance sequence
    #     sequence_id, probable_error_pos, edit_distance = search_sequences(seq)
    #     values["probable_gt"] = sequence_id
    #     values["probable_err"] = edit_distance
    #     values["probable_err_pos"] = probable_error_pos
    # else:
    #     # we found this array either in excluding primer match that means all the error is in the priemer
    #     values["probable_err"] = values['probable_err_in_primer']
    #     values["probable_err_pos"] = ""
    #     cmd = f"INSERT INTO sequences (length, r1_count, r2_count, gt_including_primer, gt_pc, gt_pl, probable_gt, " \
    #           f"probable_err, probable_err_in_primer, probable_err_pos, sequence) VALUES (:length, :r1_count, :r2_count, " \
    #           f":gt_including_primer, :gt_pc,:gt_pl, :probable_gt, :probable_err, :probable_err_in_primer, " \
    #           f":probable_err_pos, :sequence)"
    pass


def search_sequences(sequence):

    if len(sequence) < 220:
        return -1, "", -1  # consider this sequence as garbage
    cursor = connection.cursor()
    cursor.execute("SELECT id, sequence_without_primer FROM gt_sequences")
    gt_sequences = cursor.fetchall()
    # gt_id, pgt_sequence, min_edit_distance = -1, "", float('inf')
    p_get_edit_distance = partial(get_edit_distance, sequence=sequence)
    optimum_number_of_process = int(math.ceil(multiprocessing.cpu_count()))
    pool = multiprocessing.Pool(processes=optimum_number_of_process)
    edit_distances = pool.map(p_get_edit_distance, gt_sequences)
    pool.close()
    pool.join()
    # sort the edit distance
    edit_distances.sort(key=lambda x: x[1])
    # TODO: implement this
    probable_error_pos = ""
    return edit_distances[0][0], probable_error_pos, edit_distances[0][1]
    # for id_, gt_sequence in gt_sequences:
    #     edit_distance = get_edit_distance(sequence, gt_sequence)
    #     if edit_distance < min_edit_distance:
    #         min_edit_distance = edit_distance
    #         gt_id = id_
    #     if min_edit_distance == 1:
    #         break
    # Find probable error position
    # TODO: implement this
    # probable_error_pos = ""
    # return gt_id, probable_error_pos, min_edit_distance


def get_edit_distance(gt_sequence_with_id_: str, sequence) -> int:
    gt_sequence = gt_sequence_with_id_[1]
    delta = lambda i, j: 1 if sequence[i] != gt_sequence[j] else 0
    n = len(gt_sequence)
    m = len(sequence)
    dp = [[0 for i in range(n + 1)] for j in range(m + 1)]

    for i in range(1, m + 1):
        dp[i][0] = i

    for j in range(1, n + 1):
        dp[0][j] = j

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            dp[i][j] = min(
                dp[i - 1][j - 1] + delta(i - 1, j - 1),
                dp[i - 1][j] + 1,
                dp[i][j - 1] + 1,
            )  # align, delete, insert respectively

    return gt_sequence_with_id_[0], dp[m][n]


insert_wetlab_data("seqNAM01_R1_001.fastq", 1)

insert_wetlab_data("seqNAM01_R2_001.fastq", 2)

connection.close()
