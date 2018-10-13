from Bio.SeqUtils.ProtParam import ProteinAnalysis
import Bio.SeqUtils.ProtParamData

def get_amino_acids():
    return ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
    
class AminoAcid:

    # properties = ['charge','phosphorylation','average_flexibility_idx','ionic_bond','molecular_weight','hydrophobicity','typically_helix','typically_turn','typically_sheet','hydrophilicity','surface_accessibility','mutability','janin_interior_surface_energy_scale']
    properties = ['charge', 'phosphorylation', 'average_flexibility_idx',
                  'ionic_bond', 'molecular_weight', 'hydrophobicity']


    def __init__(self,amino_acid):

        self.properties = {}

        self.is_valid = True

        if amino_acid == 'R':
            self.is_valid = True
            self.amino_acid_letter = 'R'
            self.amino_acid_name = 'Arginine'
            self.properties['charge'] = 1.00
            self.properties['hydropathy'] = -4.5
            self.properties['solubility'] = -14.0
            self.properties['phosphorylation'] = 0.00
            self.properties['average_flexibility_idx'] = 0.530
            self.properties['ionic_bond'] = 1.00
            self.properties['hydrogen_bond'] = 1.00
            self.properties['hydrophilicity'] = 3.0
            self.properties['surface_accessibility'] = 1.475
            self.properties['mutability'] = 65

        elif amino_acid== 'N':
            self.is_valid = True
            self.amino_acid_letter = 'N'
            self.amino_acid_name = 'Asparagine'
            self.properties['charge'] = 0.00
            self.properties['hydropathy'] = -3.5
            self.properties['solubility'] = -28
            self.properties['phosphorylation'] = 0.00
            self.properties['average_flexibility_idx'] = 0.46
            self.properties['ionic_bond'] = 0.00
            self.properties['hydrogen_bond'] = 1.00
            self.properties['hydrophilicity'] = 0.2
            self.properties['hydrophilicity'] = 3.0
            self.properties['surface_accessibility'] = 1.296
            self.properties['mutability'] = 134

        elif amino_acid == 'D':
            self.is_valid = True
            self.amino_acid_letter = 'D'
            self.amino_acid_name = 'Aspartate'
            self.properties['charge'] = -1.00
            self.properties['hydropathy'] = -3.5
            self.properties['solubility'] = -55
            self.properties['phosphorylation'] = 0.00
            self.properties['average_flexibility_idx'] = 0.510
            self.properties['ionic_bond'] = 1.00
            self.properties['hydrogen_bond'] = 1.00
            self.properties['hydrophilicity'] = 3.0
            self.properties['surface_accessibility'] = 1.283
            self.properties['mutability'] = 106

        elif amino_acid == 'E':
            self.is_valid = True
            self.amino_acid_letter = 'E'
            self.amino_acid_name = 'Glutamate'
            self.properties['charge'] = -1.00
            self.properties['hydropathy'] = -3.5
            self.properties['solubility'] = -31
            self.properties['phosphorylation'] = 0.0
            self.properties['average_flexibility_idx'] = 0.500
            self.properties['ionic_bond'] = 1.00
            self.properties['hydrogen_bond'] = 1.00
            self.properties['hydrophilicity'] = 3.0
            self.properties['surface_accessibility'] = 1.445
            self.properties['mutability'] = 102

        elif amino_acid == 'Q':
            self.is_valid = True
            self.amino_acid_letter = 'Q'
            self.amino_acid_name = 'Glutamine'
            self.properties['charge'] = 0.00
            self.properties['hydropathy'] = -3.5
            self.properties['solubility'] = -10
            self.properties['phosphorylation'] = 0.0
            self.properties['average_flexibility_idx'] = 0.490
            self.properties['ionic_bond'] = 0.00
            self.properties['hydrogen_bond'] = 1.00
            self.properties['hydrophilicity'] = 0.2
            self.properties['surface_accessibility'] = 1.348
            self.properties['mutability'] = 93

        elif amino_acid == 'K':
            self.is_valid = True
            self.amino_acid_letter = 'K'
            self.amino_acid_name = 'Lysine'
            self.properties['charge'] = 1.00
            self.properties['hydropathy'] = -3.9
            self.properties['solubility'] = -23
            self.properties['phosphorylation'] = 0.0
            self.properties['average_flexibility_idx'] = 0.470
            self.properties['ionic_bond'] = 1.00
            self.properties['hydrogen_bond'] = 1.00
            self.properties['hydrophilicity'] = 3.0
            self.properties['surface_accessibility'] = 1.545
            self.properties['mutability'] = 56

        elif amino_acid == 'S':
            self.is_valid = True
            self.amino_acid_letter = 'S'
            self.amino_acid_name = 'Serine'
            self.properties['charge'] = 0.00
            self.properties['hydropathy'] = -0.8
            self.properties['solubility'] = -5
            self.properties['phosphorylation'] = 1.00
            self.properties['average_flexibility_idx'] = 0.510
            self.properties['ionic_bond'] = 0.00
            self.properties['hydrogen_bond'] = 1.00
            self.properties['hydrophilicity'] = 0.3
            self.properties['surface_accessibility'] = 1.115
            self.properties['mutability'] = 120

        elif amino_acid == 'T':
            self.is_valid = True
            self.amino_acid_letter = 'T'
            self.amino_acid_name = 'Threonine'
            self.properties['charge'] = 0.00
            self.properties['hydropathy'] = -0.7
            self.properties['solubility'] = 13
            self.properties['phosphorylation'] = 1.00
            self.properties['average_flexibility_idx'] = 0.440
            self.properties['ionic_bond'] = 0.00
            self.properties['hydrogen_bond'] = 1.00
            self.properties['hydrophilicity'] = -0.4
            self.properties['surface_accessibility'] = 1.184
            self.properties['mutability'] = 97

        elif amino_acid == 'C':
            self.is_valid = True
            self.amino_acid_letter = 'C'
            self.amino_acid_name = 'Cysteine'
            self.properties['charge'] = 0
            self.properties['hydropathy'] = 2.5
            self.properties['solubility'] = 49
            self.properties['phosphorylation'] = 0.0
            self.properties['average_flexibility_idx'] = 0.350
            self.properties['ionic_bond'] = 0.00
            self.properties['hydrogen_bond'] = 0.00
            self.properties['hydrophilicity'] = -1.0
            self.properties['surface_accessibility'] = 0.394
            self.properties['mutability'] = 20

        elif amino_acid == 'H':
            self.is_valid = True
            self.amino_acid_letter = 'H'
            self.amino_acid_name = 'Histidine'
            self.properties['charge'] = 1.00
            self.properties['hydropathy'] = -3.2
            self.properties['solubility'] = 8
            self.properties['phosphorylation'] = 0.0
            self.properties['average_flexibility_idx'] = 0.320
            self.properties['ionic_bond'] = 1.00
            self.properties['hydrogen_bond'] = 1.00
            self.properties['hydrophilicity'] = -0.5
            self.properties['surface_accessibility'] = 1.180
            self.properties['mutability'] = 66


        elif amino_acid == 'M':
            self.is_valid = True
            self.amino_acid_letter = 'M'
            self.amino_acid_name = 'Methionine'
            self.properties['charge'] = 0.00
            self.properties['hydropathy'] = 1.9
            self.properties['solubility'] = 74
            self.properties['phosphorylation'] = 0.0
            self.properties['average_flexibility_idx'] = 0.300
            self.properties['ionic_bond'] = 0.00
            self.properties['hydrogen_bond'] = 0.00
            self.properties['hydrophilicity'] = -1.3
            self.properties['surface_accessibility'] = 0.714
            self.properties['mutability'] = 94


        elif amino_acid == 'A':
            self.is_valid = True
            self.amino_acid_letter = 'A'
            self.amino_acid_name = 'Alanine'
            self.properties['charge'] = 0
            self.properties['hydropathy'] = 1.8
            self.properties['solubility'] = 41
            self.properties['phosphorylation'] = 0.0
            self.properties['average_flexibility_idx'] = 0.360
            self.properties['ionic_bond'] = 0.00
            self.properties['hydrogen_bond'] = 0.00
            self.properties['hydrophilicity'] = -0.5
            self.properties['surface_accessibility'] = 0.815
            self.properties['mutability'] = 100

        elif amino_acid == 'V':
            self.is_valid = True
            self.amino_acid_letter = 'V'
            self.amino_acid_name = 'Valine'
            self.properties['charge'] = 0
            self.properties['hydropathy'] = 4.2
            self.properties['solubility'] = 76
            self.properties['phosphorylation'] = 0.0
            self.properties['average_flexibility_idx'] = 0.390
            self.properties['ionic_bond'] = 0.00
            self.properties['hydrogen_bond'] = 0.00
            self.properties['hydrophilicity'] = -1.5
            self.properties['surface_accessibility'] = 0.606
            self.properties['mutability'] = 74

        elif amino_acid == 'G':
            self.is_valid = True
            self.amino_acid_letter = 'G'
            self.amino_acid_name = 'Glycine'
            self.properties['charge'] = 0
            self.properties['hydropathy'] = -0.4
            self.properties['solubility'] = 0
            self.properties['phosphorylation'] = 0.0
            self.properties['average_flexibility_idx'] = 0.540
            self.properties['ionic_bond'] = 0.00
            self.properties['hydrogen_bond'] = 0.00
            self.properties['hydrophilicity'] = 0.0
            self.properties['surface_accessibility'] = 0.714
            self.properties['mutability'] = 49

        elif amino_acid == 'I':
            self.is_valid = True
            self.amino_acid_letter = 'I'
            self.amino_acid_name = 'Isoleucine'
            self.properties['charge'] = 0
            self.properties['hydropathy'] = 4.5
            self.properties['solubility'] = 99
            self.properties['phosphorylation'] = 0.0
            self.properties['average_flexibility_idx'] = 0.460
            self.properties['ionic_bond'] = 0.00
            self.properties['hydrogen_bond'] = 0.00
            self.properties['hydrophilicity'] = -1.8
            self.properties['surface_accessibility'] = 0.603
            self.properties['mutability'] = 96

        elif amino_acid == 'L':
            self.is_valid = True
            self.amino_acid_letter = 'L'
            self.amino_acid_name = 'Leucine'
            self.properties['charge'] = 0
            self.properties['hydropathy'] = 3.8
            self.properties['solubility'] = 97
            self.properties['phosphorylation'] = 0.0
            self.properties['average_flexibility_idx'] = 0.370
            self.properties['ionic_bond'] = 0.00
            self.properties['hydrogen_bond'] = 0.00
            self.properties['hydrophilicity'] = -1.8
            self.properties['surface_accessibility'] = 0.603
            self.properties['mutability'] = 40

        elif amino_acid == 'F':
            self.is_valid = True
            self.amino_acid_letter = 'F'
            self.amino_acid_name = 'Phenylalanine'
            self.properties['charge'] = 0.0
            self.properties['hydropathy'] = 2.8
            self.properties['solubility'] = 100
            self.properties['phosphorylation'] = 0
            self.properties['average_flexibility_idx'] = 0.310
            self.properties['ionic_bond'] = 0.00
            self.properties['hydrogen_bond'] = 0.00
            self.properties['hydrophilicity'] = -2.5
            self.properties['surface_accessibility'] = 0.695
            self.properties['mutability'] = 41

        elif amino_acid =='P':
            self.is_valid = True
            self.amino_acid_letter = 'P'
            self.amino_acid_name = 'Proline'
            self.properties['charge'] = 0
            self.properties['hydropathy'] = -1.6
            self.properties['solubility'] = -46
            self.properties['phosphorylation'] = 0.0
            self.properties['average_flexibility_idx'] = 0.510
            self.properties['ionic_bond'] = 0.00
            self.properties['hydrogen_bond'] = 0.00
            self.properties['hydrophilicity'] = 0.0
            self.properties['surface_accessibility'] = 1.236
            self.properties['mutability'] = 56

        elif amino_acid == 'W':
            self.is_valid = True
            self.amino_acid_letter = 'W'
            self.amino_acid_name = 'Tryptophan'
            self.properties['charge'] = 0
            self.properties['hydropathy'] = -0.9
            self.properties['solubility'] = 97
            self.properties['phosphorylation'] = 0.0
            self.properties['average_flexibility_idx'] = 0.310
            self.properties['ionic_bond'] = 0.00
            self.properties['hydrogen_bond'] = 1.00
            self.properties['hydrophilicity'] = -3.4
            self.properties['surface_accessibility'] = 0.808
            self.properties['mutability'] = 18

        elif amino_acid == 'Y':
            self.is_valid = True
            self.amino_acid_letter = 'Y'
            self.amino_acid_name = 'Tyrosine'
            self.properties['charge'] = 0
            self.properties['hydropathy'] = -1.3
            self.properties['solubility'] = 63
            self.properties['phosphorylation'] = 1.00
            self.properties['average_flexibility_idx'] = 0.420
            self.properties['ionic_bond'] = 0.00
            self.properties['hydrogen_bond'] = 1.00
            self.properties['hydrophilicity'] = -2.3
            self.properties['surface_accessibility'] = 1.089
            self.properties['mutability'] = 41

        else:
            self.is_valid = False
            print("Invalid Amino Acid "+amino_acid)

        if (self.is_valid):
            self.properties['molecular_weight'] = ProteinAnalysis(str(amino_acid)).molecular_weight()
            self.properties['hydrophobicity'] = ProteinAnalysis(str(amino_acid)).gravy()
            secondary_struct = ProteinAnalysis(str(amino_acid)).secondary_structure_fraction()
            self.properties['typically_helix'] = secondary_struct[0]
            self.properties['typically_turn'] = secondary_struct[1]
            self.properties['typically_sheet'] = secondary_struct[2]
            self.properties['janin_interior_surface_energy_scale'] = Bio.SeqUtils.ProtParamData.ja[str(amino_acid)]




def get_amino_acid_characteristics_features(sequence_matrix):

    feature_descriptions = []

    properties = AminoAcid.properties

    for feature_type in properties:
        for position_index in range(len(sequence_matrix[0])):
            feature_descriptions.append(feature_type + ' of position ' + str(position_index + 1))

    num_examples = len(sequence_matrix)
    num_properties = len(properties)
    num_amino_acids = len(sequence_matrix[0])

    feature_matrices = {}

    amino_acids = {}

    for amino_acid_letter in get_amino_acids():
        amino_acids[amino_acid_letter] = AminoAcid(amino_acid_letter)

    for property_index in range(0, num_properties):
        feature_matrices[property_index] = numpy.zeros([num_examples, num_amino_acids], dtype=numpy.float32)

        min_value = list(amino_acids.values())[0].properties[properties[property_index]]
        max_value = list(amino_acids.values())[0].properties[properties[property_index]]

        for amino_acid_letter, amino_acid in amino_acids.items():
            property_value = amino_acid.properties[properties[property_index]]
            min_value = min(property_value, min_value)
            max_value = max(property_value, max_value)

        for amino_acid_letter, amino_acid in amino_acids.items():
            property_value = amino_acid.properties[properties[property_index]]
            amino_acid.properties[properties[property_index]] = (property_value - min_value) * 1.0 / (max_value - min_value)

    for sequence_index in range(0, num_examples):

        for amino_acid_index in range(0, num_amino_acids):

            for property_index in range(0, num_properties):

                feature_matrices[property_index][sequence_index][amino_acid_index] = amino_acids[sequence_matrix[sequence_index][amino_acid_index]].properties[properties[property_index]]

    feature_matrix = numpy.zeros([num_examples, 0], dtype=numpy.float32)

    for property_index in range(0, num_properties):

        feature_matrix = numpy.concatenate((feature_matrix, feature_matrices[property_index]), axis=1)

    return feature_descriptions, feature_matrix, properties
