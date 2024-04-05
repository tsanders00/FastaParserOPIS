package ProjektOPIS;

import java.util.Dictionary;
import java.util.Enumeration;
import java.util.Hashtable;

import static java.lang.Math.round;

public class RNA extends AbstractFastaObjekt {
    private int length;
    private Dictionary<Character, Integer> composition;
    private String sequence;
    private int ID;
    private sequencetype sequencetype;
    private double gcContent;
    private double molecularWeight;

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
        dict.put('U', 0);
        dict.put('G', 0);

        for (int i = 0; i < sequence.length(); i++) {
            char b = sequence.charAt(i);
            for (Enumeration<Character> enm = dict.keys(); enm.hasMoreElements();) {
                char nextElement = enm.nextElement();
                if (b == nextElement) {
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
     * @return molecular weight of RNA seq
     */
    public double molecularWeight() {
        double A_MG = 329.2;
        double C_MG = 305.2;
        double G_MG = 345.2;
        double U_MG = 306.2;
        double molekulargewicht;
        molekulargewicht = ((composition.get('A')) * A_MG) + ((composition.get('C')) * C_MG) +
                ((composition.get('G')) * G_MG) + ((composition.get('U')) * U_MG) + 159;
        this.molecularWeight = round(molekulargewicht);
        return round(molekulargewicht);
    }
    /**
     * method calculates the percentage of how often G and C appear within the sequence
     *
     * @return gc percentage of RNA seq
     */
    public double gcContent() {
        double gcGehalt;
        gcGehalt = (double) ((composition.get('G')) + (composition.get('C'))) / (double) this.length;
        this.gcContent = gcGehalt * 100;
        return gcGehalt * 100;
    }
    /**
     * method to print basic information for a certain fasta object
     *
     * @return basic information about fasta object
     */
    public String getInformation() {
        return ("This Fasta object of type " + this.sequencetype + " with the sequence " + this.sequence
                + " and the ID " + this.ID + " has a length of " + this.length + "bp. Further details: GC-content: "
                + this.gcContent + "%, molecular weight in Dalton: " + this.molecularWeight);
    }

    @Override
    public void meltingTemp() {
        //stays empty
    }

    @Override
    public String nucleotoamino() {
        return null;
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