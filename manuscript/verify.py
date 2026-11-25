import itertools
import json
from fractions import Fraction


def require(condition, message):
    if not condition:
        raise ValueError(message)


def labelled_distances(marks):
    return sorted((marks[right] - marks[left], left, right)
                  for left in range(len(marks))
                  for right in range(left + 1, len(marks)))


def perfect_difference_set(marks, modulus):
    residues = set()
    for left in marks:
        for right in marks:
            if left == right:
                continue
            difference = (right - left) % modulus
            if difference == 0 or difference in residues:
                return False
            residues.add(difference)
    return len(residues) == modulus - 1


def interval_row(left, right, dimension):
    return tuple(int(left <= index < right) for index in range(dimension))


def sparse_certificate(marks):
    dimension = len(marks) - 1
    gaps = tuple(right - left for left, right in zip(marks, marks[1:]))
    distances = labelled_distances(marks)
    require(len({item[0] for item in distances}) == len(distances), 'Repeated distance')
    by_distance = {distance: (left, right) for distance, left, right in distances}
    require(all(distance in by_distance for distance in range(1, max(gaps) + 1)),
            'Incomplete distance prefix')
    predecessors = [interval_row(*by_distance[gap - 1], dimension) if gap > 1
                    else tuple(0 for column in range(dimension)) for gap in gaps]
    rows = [tuple(int(index == column) - entry for column, entry in enumerate(predecessor))
            for index, predecessor in enumerate(predecessors)]
    order = sorted(range(dimension), key=lambda index: gaps[index])
    for index in range(dimension):
        require(all(not entry or gaps[column] < gaps[index]
                    for column, entry in enumerate(predecessors[index])),
                'Sparse matrix is not unit triangular in gap order')
        require(sum(entry * gap for entry, gap in zip(rows[index], gaps)) == 1,
                'Sparse row is not tight')
    weights = [0] * dimension
    for index in reversed(order):
        weights[index] = 1 + sum(weights[other] * predecessors[other][index]
                                 for other in range(dimension))
    require(all(weight >= 1 for weight in weights), 'Sparse weight is not positive')
    require(all(sum(weights[index] * rows[index][column] for index in range(dimension)) == 1
                for column in range(dimension)), 'Sparse dual objective differs')
    require(sum(weights) == marks[-1], 'Sparse dual value differs')
    return {'rows': [list(row) for row in rows], 'weights': weights,
            'gap_order': order, 'determinant_in_gap_order': 1, 'bound': sum(weights)}


def certificate(marks, cut_gap):
    dimension = len(marks) - 1
    gaps = tuple(right - left for left, right in zip(marks, marks[1:]))
    maximum = max(gaps)
    distances = labelled_distances(marks)
    require(len({item[0] for item in distances}) == len(distances), 'Repeated distance')
    require(maximum <= cut_gap, 'Cut is not at a largest gap')
    by_distance = {distance: (left, right) for distance, left, right in distances}
    require(all(distance in by_distance for distance in range(1, maximum + 1)),
            'Incomplete prefix of ordinary distances')
    vectors = [tuple(0 for index in range(dimension))]
    for distance in range(1, maximum + 1):
        vectors.append(interval_row(*by_distance[distance], dimension))
    rows = [tuple(current - previous for current, previous in zip(vectors[distance],
                                                                 vectors[distance - 1]))
            for distance in range(1, maximum + 1)]
    weights = [sum(gap >= distance for gap in gaps) for distance in range(1, maximum + 1)]
    objective = tuple(sum(weight * row[index] for weight, row in zip(weights, rows))
                      for index in range(dimension))
    require(objective == tuple(1 for index in range(dimension)), 'Wrong dual objective')
    require(sum(weights) == marks[-1], 'Wrong dual value')
    require(all(weight > 0 for weight in weights), 'Nonpositive weight')
    for row in rows:
        require(sum(entry * gap for entry, gap in zip(row, gaps)) == 1, 'Row is not tight')
    for index, gap in enumerate(gaps):
        cumulative = tuple(sum(row[column] for row in rows[:gap])
                           for column in range(dimension))
        require(cumulative == tuple(int(index == column) for column in range(dimension)),
                'Tight rows do not recover the gap coordinates')
    return {'marks': list(marks), 'cut_gap': cut_gap, 'gaps': list(gaps),
            'rows': [list(row) for row in rows], 'weights': weights,
            'objective': list(objective), 'bound': sum(weights),
            'sparse_certificate': sparse_certificate(marks)}


def test_perfect_sets():
    counts = []
    examples = []
    for cardinality in range(2, 7):
        modulus = cardinality * (cardinality - 1) + 1
        perfect_sets = 0
        cuts = 0
        for tail in itertools.combinations(range(1, modulus), cardinality - 1):
            marks = (0,) + tail
            if not perfect_difference_set(marks, modulus):
                continue
            perfect_sets += 1
            cyclic_gaps = [right - left for left, right in zip(marks, marks[1:])]
            cyclic_gaps.append(modulus - marks[-1])
            require(len(set(cyclic_gaps)) == cardinality, 'Cyclic gaps are not distinct')
            largest = max(cyclic_gaps)
            for gap_index, gap in enumerate(cyclic_gaps):
                if gap != largest:
                    continue
                start = marks[(gap_index + 1) % cardinality]
                normalized = tuple(sorted((mark - start) % modulus for mark in marks))
                require(max(right - left for left, right in zip(normalized, normalized[1:])) < largest,
                        'Retained gap is not strictly smaller than the removed gap')
                require(normalized[-1] + largest == modulus, 'Wrong cut normalization')
                current = certificate(normalized, largest)
                cuts += 1
                if not any(example['modulus'] == modulus for example in examples):
                    examples.append({'modulus': modulus, **current})
        require(perfect_sets > 0, 'Empty perfect-set family')
        counts.append({'cardinality': cardinality, 'modulus': modulus,
                       'sets_containing_zero': perfect_sets, 'largest_gap_cuts_checked': cuts})
    return counts, examples


def test_nonlargest_cut():
    reference = (0, 4, 24, 27, 67, 77, 86, 88, 103, 119, 125, 132)
    doubled_candidate = (0, 8, 44, 50, 123, 141, 157, 161, 188, 217, 228, 242)
    candidate = tuple(Fraction(mark, 2) for mark in doubled_candidate)
    reference_order = labelled_distances(reference)
    candidate_order = labelled_distances(candidate)
    require(perfect_difference_set(reference, 133), 'Not a cyclic perfect difference set')
    require([(left, right) for distance, left, right in reference_order] ==
            [(left, right) for distance, left, right in candidate_order], 'Different chamber')
    values = [distance for distance, left, right in candidate_order]
    require(values[0] >= 1, 'Smallest distance is below one')
    require(all(right - left >= 1 for left, right in zip(values, values[1:])),
            'Distances are not one-separated')
    require(candidate[-1] < reference[-1], 'No shortening')
    return {'modulus': 133, 'reference_marks': list(reference),
            'candidate_marks': [str(mark) for mark in candidate],
            'cut_gap': 1, 'largest_cyclic_gap': 40,
            'reference_length': 132, 'candidate_length': '121',
            'source': 'Notik, arXiv:2611.22656v1, Proposition 2, pp. 5–6'}


def check_shortening(reference, candidate):
    reference_order = labelled_distances(reference)
    require(len({row[0] for row in reference_order}) == len(reference_order),
            'Reference is not Golomb')
    values = [candidate[right] - candidate[left] for distance, left, right in reference_order]
    require(values[0] >= 1 and all(right - left >= 1 for left, right in zip(values, values[1:])),
            'Candidate is outside the unit-separated reference chamber')
    require(candidate[-1] < reference[-1], 'Candidate is not shorter')
    return {'reference_marks': list(reference), 'candidate_marks': [str(mark) for mark in candidate],
            'reference_ordered_distances': [row[0] for row in reference_order],
            'candidate_distances_in_reference_order': [str(value) for value in values],
            'minimum_distance_increment': str(min(right - left for left, right in zip(values, values[1:])))}


def test_bose_ruzsa():
    powers = []
    value = (1, 0)
    for exponent in range(48):
        require(value not in powers, 'Bose generator repeats before order 48')
        powers.append(value)
        value = ((value[0] + 3 * value[1]) % 7, (value[0] + value[1]) % 7)
    require(value == (1, 0), 'Bose generator does not close at 48')
    require(powers[16] == (4, 0) and powers[24] == (6, 0), 'Bose order witnesses differ')
    require(3 not in {number * number % 7 for number in range(7)}, 'Bose polynomial reducible')
    exponents = [powers.index(((1 + parameter) % 7, 1)) for parameter in range(7)]
    require(exponents == [1, 18, 21, 27, 6, 31, 20], 'Bose line exponents differ')
    bose_residues = sorted(5 * exponent % 48 for exponent in exponents)
    alternate_generator = powers[11]
    require(len({powers[11 * exponent % 48] for exponent in range(48)}) == 48,
            'Alternative Bose generator is not primitive')
    alternate_exponents = sorted(next(exponent for exponent in range(48)
                                     if powers[11 * exponent % 48] ==
                                     ((alternate_generator[0] + parameter) % 7,
                                      alternate_generator[1]))
                                for parameter in range(7))
    require(alternate_exponents == [1, 3, 22, 31, 34, 44, 45],
            'Alternative Bose line exponents differ')
    require(sorted((8 + exponent) % 48 for exponent in alternate_exponents) == bose_residues,
            'Bose multiplier-free translate identity differs')
    bose_marks = tuple(sorted((residue - 30) % 48 for residue in bose_residues))
    require(bose_marks == (0, 9, 12, 22, 23, 27, 29), 'Bose cut differs')
    require(max(right - left for left, right in zip(bose_marks, bose_marks[1:])) <= 48 - bose_marks[-1],
            'Bose cut is not largest')
    bose = check_shortening(bose_marks, (0, Fraction(17, 2), Fraction(23, 2), 21, 22, 26, 28))
    bose.update(modulus=48, cut_gap=19, exponents=exponents, residues=bose_residues,
                missing_prefix_distance=8, generator_order=48)
    require(8 not in {row[0] for row in labelled_distances(bose_marks)}, 'Bose missing distance occurs')
    ruzsa_residues = sorted(next(residue for residue in range(20)
                                 if residue % 4 == index and residue % 5 == pow(2, index, 5))
                            for index in range(4))
    require(len({pow(2, index, 5) for index in range(4)}) == 4, 'Ruzsa root not primitive')
    ruzsa_marks = tuple(sorted((residue - 14) % 20 for residue in ruzsa_residues))
    require(ruzsa_marks == (0, 2, 3, 9), 'Ruzsa cut differs')
    require(max(right - left for left, right in zip(ruzsa_marks, ruzsa_marks[1:])) <= 20 - ruzsa_marks[-1],
            'Ruzsa cut is not largest')
    ruzsa = check_shortening(ruzsa_marks, (0, 2, 3, 8))
    ruzsa.update(modulus=20, cut_gap=11, residues=ruzsa_residues, missing_prefix_distance=4)
    require(4 not in {row[0] for row in labelled_distances(ruzsa_marks)}, 'Ruzsa missing distance occurs')
    return {'Bose': bose, 'Ruzsa': ruzsa, 'source_question': 'Notik, arXiv:2611.22656v1, Section 10 p.8',
            'conclusion': 'Answered: no, for each universal family question.'}


def test_general_prefix():
    counts = []
    for cardinality in range(2, 6):
        checked = 0
        for tail in itertools.combinations(range(1, 19), cardinality - 1):
            marks = (0,) + tail
            distances = labelled_distances(marks)
            values = {row[0] for row in distances}
            if len(values) != len(distances):
                continue
            gaps = tuple(right - left for left, right in zip(marks, marks[1:]))
            if not set(range(1, max(gaps) + 1)) <= values:
                continue
            sparse_certificate(marks)
            candidate = tuple(Fraction(5, 4) * mark for mark in marks)
            excesses = [candidate[index + 1] - candidate[index] - gap
                        for index, gap in enumerate(gaps)]
            epsilon = candidate[-1] - marks[-1]
            require(all(excess >= 0 for excess in excesses) and sum(excesses) == epsilon,
                    'Gap stability identity differs')
            require(all(0 <= candidate[index] - mark <= epsilon for index, mark in enumerate(marks)),
                    'Coordinate stability bound differs')
            checked += 1
        counts.append({'cardinality': cardinality, 'maximum_length': 18, 'prefix_rulers_checked': checked})
    return counts


def main():
    counts, examples = test_perfect_sets()
    require(not perfect_difference_set((0, 1, 2), 7), 'Invalid difference set accepted')
    try:
        certificate((0, 1, 4, 6), 1)
    except ValueError:
        pass
    else:
        raise ValueError('Non-largest cut accepted')
    result = {'arithmetic': 'exact integers and fractions',
              'finite_scope': 'All cyclic perfect difference sets containing zero, cardinalities 2–6, '
                              'moduli 3,7,13,21,31; every largest-gap cut of each set.',
              'perfect_set_checks': counts, 'certificate_examples': examples,
              'nonlargest_cut': test_nonlargest_cut(),
              'general_complete_prefix_scope': test_general_prefix(),
              'Bose_and_Ruzsa_questions': test_bose_ruzsa(),
              'negative_controls': 'Repeated cyclic differences and a non-largest cut are rejected.',
              'status': 'All finite checks passed.'}
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == '__main__':
    main()
