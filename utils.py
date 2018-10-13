def convert_amino_acid_ids_to_letters(number_sequences):
        
    amino_acids = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
    
    amino_acid_sequences = []
    
    for number_sequence in number_sequences:

        amino_acid_sequence = "".join([amino_acids[i - 1] for i in number_sequence])

        amino_acid_sequences.append(amino_acid_sequence)
    
    return amino_acid_sequences

def convert_amino_acid_sequence_to_ids(amino_acid_sequence):
        
    amino_acids = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
    
    amino_acid_dict = {}
    
    for amino_acid_index, amino_acid in enumerate(amino_acids):
        amino_acid_dict[amino_acid] = amino_acid_index + 1
        
    number_sequence = []
    
    for amino_acid in amino_acid_sequence:
        
        number_sequence.append(amino_acid_dict[amino_acid])
    
    return number_sequence
