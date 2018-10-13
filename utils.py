def convert_amino_acid_ids_to_letters(number_sequences):
    
    amino_acid_sequences = []
    
    for number_sequence in number_sequences:
        
        amino_acids = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

        amino_acid_sequence = [amino_acids[i - 1] for i in number_sequence]

        amino_acid_sequences.append(amino_acid_sequence)
    
    return amino_acid_sequences
