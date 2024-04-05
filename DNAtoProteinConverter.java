package ProjektOPIS;

public class DNAtoProteinConverter {
    /**
     * method which translates DNA to protein sequence
     * @param dnaSequence method receives the full dna sequence
     * @return is the fully translated protein sequence
     */
    public static String translateDNAToProtein(String dnaSequence) {
        // eval if length == 3, otherwise error
        if (dnaSequence.length() % 3 != 0) {
            throw new IllegalArgumentException("The DNA sequence has to be a multiple of 3!");
        }
        // take 3 nucleotides and translates them into amino acid, appends to new protein sequence
        StringBuilder proteinSequence = new StringBuilder();
        for (int i = 0; i < dnaSequence.length(); i += 3) {
            String codon = dnaSequence.substring(i, i + 3);
            String aminoAcid = translateCodonToAminoAcid(codon);
            proteinSequence.append(aminoAcid);
        }

        return proteinSequence.toString();
    }

    /**
     * method which just resembles the code sun and is used within 'translateDNAtoProtein'
     * @param codon gets given a 3 letter codon
     * @return returns the corresponding amino acid
     */
    private static String translateCodonToAminoAcid(String codon) {
        // code sun used to translate DNA to Protein
        return switch (codon) {
            case "AAA" -> "K";
            case "AAC" -> "N";
            case "AAG" -> "K";
            case "AAU" -> "N";
            case "ACA" -> "T";
            case "ACC" -> "T";
            case "ACG" -> "T";
            case "ACU" -> "T";
            case "AGA" -> "R";
            case "AGC" -> "S";
            case "AGG" -> "R";
            case "AGU" -> "S";
            case "AUA" -> "I";
            case "AUC" -> "I";
            case "AUG" -> "M"; // start codon
            case "AUU" -> "I";
            case "CAA" -> "Q";
            case "CAC", "CAU" -> "H";
            case "CAG" -> "Q";
            case "CCG", "CCC", "CCA", "CCU" -> "P";
            case "CGA", "CGC", "CGG", "CGU" -> "R";
            case "CUA", "CUC", "CUG", "CUU", "UUA", "UUG" -> "L";
            case "GAA" -> "E";
            case "GAC", "GAU" -> "D";
            case "GAG" -> "E";
            case "GCA", "GCC", "GCG", "GCU" -> "A";
            case "GGA", "GGC", "GGG", "GGU" -> "G";
            case "GUA", "GUC", "GUG", "GUU" -> "V";
            case "UAA", "UAG", "UGA" -> "*"; // stopp-codons
            case "UAC", "UAU" -> "Y";
            case "UCA", "UCC", "UCG", "UCU" -> "S";
            case "UGC", "UGU" -> "C";
            case "UGG" -> "W";
            case "UUC", "UUU" -> "F";
            default -> "X"; // Unbekanntes Codon
        };
    }
}

