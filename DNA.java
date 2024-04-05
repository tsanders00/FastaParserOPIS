package ProjektOPIS;

import java.util.Dictionary;
import java.util.Enumeration;
import java.util.Hashtable;

import static java.lang.Math.round;

public class DNA extends AbstractFastaObjekt {

    private int length;
    private Dictionary<Character, Integer> composition;
    private String sequence;
    private String proteinSeq;
    private int ID;
    private sequencetype sequencetype;
    private double gcContent;
    private double molecularWeight;
    private double meltingTemp;

    /**
     method to add sequence to an object
     @param sequence method needs a sequence (String)
     */
    public void addSeq(String sequence) {
        this.sequence = sequence;
    }
    /**
     method to set the sequence_type
     @param typ method needs a sequence type (Enum)
     */
    public void setSequenceType(sequencetype typ) {this.sequencetype = typ;}

    /**
     method to set an ID
     @param ID the ID is the corresponding number in the array
     */
    public void addID(int ID) {this.ID = ID;}

    /**
     * method to determine the amino acids within the sequence
     * the numbers are stored in a dictionary
     */
    public void composition() {
        Dictionary<Character, Integer> dict = new Hashtable<>();
        dict.put('A', 0);
        dict.put('C', 0);
        dict.put('T', 0);
        dict.put('G', 0);

        for (int i = 0; i < sequence.length(); i++) {
            char b = sequence.charAt(i);
            for (Enumeration<Character> enm = dict.keys(); enm.hasMoreElements();) {
                char nextElement = enm.nextElement();
                if ( b == nextElement) {
                    int update = dict.get(nextElement);
                    dict.put(b, update+1);
                }
            }
        }
        this.composition = dict;
        this.length = sequence.length();
    }

    /**
     * method calculates the molecular weight of the given DNA sequence
     *
     * @return molecular weight of DNA seq
     */
    public double molecularWeight() {
        double A_MG = 313.21;
        double C_MG = 289.18;
        double G_MG = 329.21;
        double T_MG = 304.20;
        double molekulargewicht;
        molekulargewicht = ((composition.get('A')) * A_MG) + ((composition.get('C')) * C_MG)
                + ((composition.get('G')) * G_MG) + ((composition.get('T')) * T_MG) - 61.96;
        this.molecularWeight = round(molekulargewicht);
        return round(molekulargewicht);
    }

    /**
     * method calculates the percentage of how often G and C appear within the sequence
     *
     * @return gc percentage of DNA seq
     */
    public double gcContent() {
        double gcGehalt;
        gcGehalt = (double) ((composition.get('G')) + (composition.get('C'))) / (double) this.length;
        this.gcContent = gcGehalt * 100;
        return gcGehalt * 100;
    }

    /**
     * method calculates the melting temperature for a given DNA sequence
     * the calculation depends on the length of the sequence
     */
    public void meltingTemp() {
        if (this.length <14) {
            this.meltingTemp = round((((composition.get('A')) + (composition.get('T'))) * 2)
                    + ((((composition.get('G'))) + (composition.get('C'))) *4));
        }
        else {
            this.meltingTemp = round(64.9 + (41*((composition.get('G'))) + (composition.get('C')) - 16.4)
                    / ((composition.get('A')) + (composition.get('G')) +
                    (composition.get('C')) + (composition.get('T'))));
        }
    }

    /**
     * method to print basic information for a certain fasta object
     *
     * @return basic information about fasta object
     */
    public String getInformation() {
        return ("This Fasta object of type " + this.sequencetype + " with the sequence " + this.sequence + " and the ID " + this.ID +
                " has a length of " + this.length + "bp. Further details: GC-content: " + this.gcContent +
                "%, molecular weight in Dalton: " + this.molecularWeight + " , melting temperature (°C): " + this.meltingTemp);
    }

    /**
     * method which translates DNA into protein sequence
     * @return is a protein sequence
     */
    public String nucleotoamino() {
        String proteinsequence = DNAtoProteinConverter.translateDNAToProtein(this.sequence);
        this.proteinSeq = proteinsequence;
        return proteinsequence;
    }

    @Override
    public double netCharge(double ph) {
    return 0;
    }

    @Override
    public double isoelectricPoint() {
        return 0;
    }
}