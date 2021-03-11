import sqlite3
from collections import defaultdict, Counter
import math


def round_nearest_10(number):
    if number % 10 < 5:
        # round down
        return int(number/10) * 10
    else:
        # round up
        return int((number + 10) / 10) * 10


def analyze_seq_length(length_counter):
    # Right now this function uses read which have the length of the sequences.
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
    for length, count in length_counter:
        length_rounded = int(math.ceil(int(length)/10.0)) * 10
        if length_rounded == 250:
            counter[str(length)] += count
        elif length_rounded < 210:
            counter[counter_map["less"]] += count
        elif length_rounded > 290:
            counter[counter_map["greater"]] += count
        else:
            counter[counter_map[str(length_rounded)]] += count
    return dict(counter)


connection = sqlite3.connect("sequences.db")
cursor = connection.cursor()


"""Sequences that we were never able to recover"""
cursor.execute("SELECT DISTINCT gt_including_primer FROM sequences WHERE gt_including_primer > 0")
recovered_id = set(list(zip(*cursor.fetchall()))[0])
not_recovered_id = set(range(1, 605)).difference(recovered_id)
cursor.execute(f"select sequence_with_primer from gt_sequences where id IN {tuple(not_recovered_id)}")
print(f"unrecovred_sequences==>{list(zip(*cursor.fetchall()))[0]}")

# Analyze length of the sequences
# length of sequences for r1 count
cursor.execute("Select length, sum(r1_count) from sequences group by length")
lengths = cursor.fetchall()
print(f"r1_length_counter==> {analyze_seq_length(lengths)}")
cursor.execute("Select length, sum(r2_count) from sequences group by length")
lengths = cursor.fetchall()
print(f"r2_length_counter==> {analyze_seq_length(lengths)}")
# last id is the total number of unique sequence
cursor.execute("select id from sequences order by id desc limit 1")
total_sequences = cursor.fetchall()[0][0]
print(f"total unique==> {total_sequences}")
# sequence in read 1
cursor.execute("select count(id) from sequences where r1_count > 0")
total_sequences_read_1 = cursor.fetchall()[0][0]
print(f"r1 total unique==> {total_sequences_read_1}")
# sequence in read 2
cursor.execute("select count(id) from sequences where r2_count > 0")
total_sequences_read_2 = cursor.fetchall()[0][0]
print(f"r2 total unique==> {total_sequences_read_2}")

"""
Analyze correctness
"""
# total correctness including primer
cursor.execute("select r1_count, r2_count from sequences where gt_including_primer > 0")
total_correct_sequences = cursor.fetchall()
print(f"total unique correct==> {len(total_correct_sequences)}")
r1_correct, r2_correct = 0, 0
r1_unique_correct, r2_unique_correct = 0, 0
for r1, r2 in total_correct_sequences:
    if r1:
        r1_correct += r1
        r1_unique_correct += 1
    if r2:
        r2_correct += r2
        r2_unique_correct += 1
print(f"total  correct==> {r1_correct + r2_correct}")
# read one
print(f"r1 unique correct==> {r1_unique_correct}")
print(f"r2 unique correct==> {r2_unique_correct}")
print(f"r1 total correct==> {r1_correct}")
print(f"r2 total correct==> {r2_correct}")

# excluding primer based on primer content
cursor.execute("select r1_count, r2_count from sequences where gt_pc > 0")
total_correct_sequences = cursor.fetchall()
print(f"total unique correct(excluding primer based on content)==> {len(total_correct_sequences)}")
r1_correct, r2_correct = 0, 0
r1_unique_correct, r2_unique_correct = 0, 0
for r1, r2 in total_correct_sequences:
    if r1:
        r1_correct += r1
        r1_unique_correct += 1
    if r2:
        r2_correct += r2
        r2_unique_correct += 1
print(f"total correct(excluding primer based on content)==> {r1_correct + r2_correct}")
# read one
print(f"r1 unique correct(excluding primer based on content)==> {r1_unique_correct}")
print(f"r2 unique correct(excluding primer based on content)==> {r2_unique_correct}")
print(f"r1 total correct(excluding primer based on content)==> {r1_correct}")
print(f"r2 total correct(excluding primer based on content)==> {r2_correct}")

# excluding primer based on primer content
cursor.execute("select r1_count, r2_count from sequences where gt_pl > 0")
total_correct_sequences = cursor.fetchall()
print(
    f"total unique correct sequences in both reads(excluding primer based on length)==> {len(total_correct_sequences)}")
r1_correct, r2_correct = 0, 0
r1_unique_correct, r2_unique_correct = 0, 0
for r1, r2 in total_correct_sequences:
    if r1:
        r1_correct += r1
        r1_unique_correct += 1
    if r2:
        r2_correct += r2
        r2_unique_correct += 1
print(f"total  correct sequences in both reads(excluding primer based on length)==> {r1_correct + r2_correct}")
# read one
print(f"r1 unique correct(excluding primer based on length)==> {r1_unique_correct}")
print(f"r2 unique correct(excluding primer based on length)==> {r2_unique_correct}")
print(f"r1 total correct(excluding primer based on length)==> {r1_correct}")
print(f"r2 total correct(excluding primer based on length)==> {r2_correct}")

"""
Error distribution
"""
cursor.execute("select r1_count, r2_count, probable_err, probable_err_in_primer from sequences WHERE probable_gt > 0")
result = cursor.fetchall()

unique_error_distribution = {}
r1_unique_err_distribution = {}
r2_unique_err_distribution = {}

unique_primer_error_distribution = {}
r1_unique_primer_err_distribution = {}
r2_unique_primer_err_distribution = {}

total_error_distribution = {}
r1_total_err_distribution = {}
r2_total_err_distribution = {}

total_primer_error_distribution = {}
r1_total_primer_err_distribution = {}
r2_total_primer_err_distribution = {}


for single_result in result:
    if single_result[0]:  # r1_count
        r1_unique_err_distribution[single_result[2]] = r1_unique_err_distribution.get(single_result[2], 0) + 1
        r1_total_err_distribution[single_result[2]] = r1_total_err_distribution.get(single_result[2], 0) + \
                                                      single_result[0]
        total_error_distribution[single_result[2]] = total_error_distribution.get(single_result[2], 0) + \
                                                     single_result[0]

        r1_unique_primer_err_distribution[single_result[3]] = r1_unique_primer_err_distribution.get(single_result[3],
                                                                                                    0) + 1
        r1_total_primer_err_distribution[single_result[3]] = r1_total_primer_err_distribution.get(single_result[3], 0) + \
                                                             single_result[0]
        total_primer_error_distribution[single_result[3]] = total_primer_error_distribution.get(single_result[3], 0) + \
                                                            single_result[0]

    if single_result[1]:  # r2_count
        r2_unique_err_distribution[single_result[2]] = r2_unique_err_distribution.get(single_result[2], 0) + 1
        r2_total_err_distribution[single_result[2]] = r2_total_err_distribution.get(single_result[2], 0) + \
                                                      single_result[1]
        total_error_distribution[single_result[2]] = total_error_distribution.get(single_result[2], 0) + \
                                                     single_result[1]

        r2_unique_primer_err_distribution[single_result[3]] = r2_unique_primer_err_distribution.get(single_result[3],
                                                                                                    0) + 1
        r2_total_primer_err_distribution[single_result[3]] = r2_total_primer_err_distribution.get(single_result[3], 0) + \
                                                             single_result[0]
        total_primer_error_distribution[single_result[3]] = total_primer_error_distribution.get(single_result[3], 0) + \
                                                            single_result[0]

    unique_error_distribution[single_result[2]] = unique_error_distribution.get(single_result[2], 0) + 1
    unique_primer_error_distribution[single_result[3]] = unique_error_distribution.get(single_result[3], 0) + 1
del result
# error distribution in total sequences
print(f"unique error distribution==>{unique_error_distribution}")
print(f"r1_unique_err_distribution==>{r1_unique_err_distribution}")
print(f"r2_unique_err_distribution==>{r2_unique_err_distribution}")
print(f"unique_primer_error_distribution==>{unique_primer_error_distribution}")
print(f"r1_unique_primer_err_distribution==>{r1_unique_primer_err_distribution}")
print(f"r2_unique_primer_err_distribution==>{r2_unique_primer_err_distribution}")
print(f"total_error_distribution==>{total_error_distribution}")
print(f"r1_total_err_distribution==>{r1_total_err_distribution}")
print(f"r2_total_err_distribution==>{r2_total_err_distribution}")
print(f"total_primer_error_distribution==>{total_primer_error_distribution}")
print(f"r1_total_primer_err_distribution==>{r1_total_primer_err_distribution}")
print(f"r2_total_primer_err_distribution==>{r2_total_primer_err_distribution}")

"""sequence frequency analysis"""
# cursor.execute("Select r1_count, count(*) from sequences where r1_count > 0 group by r1_count")
cursor.execute("SELECT id, r1_count FROM sequences WHERE r1_count > 0")
r1_frequency = dict(cursor.fetchall())
total_sequences_in_r1 = sum(r1_frequency.values())
r1_frequency = Counter(r1_frequency.values())
print(f"r1 total sequences==> {total_sequences_in_r1}")
# r1_frequency = {k: round(v / total_sequences_in_r1, 2) * 100 for k, v in r1_frequency.items()}
print(f"r1 sequence frequency==> {r1_frequency}")

# cursor.execute("Select r2_count, count(*) from sequences group by r2_count")
cursor.execute("SELECT id, r2_count FROM sequences WHERE r2_count > 0")
r2_frequency = dict(cursor.fetchall())
total_sequences_in_r2 = sum(r2_frequency.values())
r2_frequency = Counter(r2_frequency.values())
print(f"r2 total sequences==> {total_sequences_in_r2}")
# r2_frequency = {k: round(v / total_sequences_in_r2, 2) * 100 for k, v in r2_frequency.items()}
print(f"r2 sequence frequency==> {r2_frequency}")
print(f"total sequence==>{total_sequences_in_r1 + total_sequences_in_r2}")
"""frequency of our ground truth data"""
cursor.execute("select gt_including_primer, r1_count, r2_count from sequences where gt_including_primer > 0")
result = cursor.fetchall()
total_gt_frequency, r1_gt_frequency, r2_gt_frequency = {}, {}, {}
for single_result in result:
    total_gt_frequency[single_result[0]] = single_result[1] + single_result[2]
    r1_gt_frequency[single_result[0]] = single_result[1]
    r2_gt_frequency[single_result[0]] = single_result[2]
print(f"total_gt_frequency==>{total_gt_frequency}")
print(f"r1_gt_frequency==>{r1_gt_frequency}")
print(f"r2_gt_frequency==>{r2_gt_frequency}")


connection.close()
