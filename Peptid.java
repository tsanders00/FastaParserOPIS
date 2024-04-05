package ProjektOPIS;

import java.util.Dictionary;
import java.util.Enumeration;
import java.util.Hashtable;

public class Peptid extends AbstractFastaObjekt {
    private Dictionary<Character, Integer> composition;
    private String sequence;
    private sequencetype sequencetype;
    private int ID;
    private int length;
    private double netCharge;
    private double isoelectricPoint;

    public Peptid() {
    }
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
        dict.put('D', 0);
        dict.put('E', 0);
        dict.put('F', 0);
        dict.put('G', 0);
        dict.put('H', 0);
        dict.put('I', 0);
        dict.put('K', 0);
        dict.put('L', 0);
        dict.put('M', 0);
        dict.put('N', 0);
        dict.put('O', 0);
        dict.put('P', 0);
        dict.put('Q', 0);
        dict.put('R', 0);
        dict.put('S', 0);
        dict.put('T', 0);
        dict.put('V', 0);
        dict.put('W', 0);
        dict.put('Y', 0);
        for (int i = 0; i < sequence.length(); i++) {
            char b = sequence.charAt(i);
            //System.out.println(b);
            for (Enumeration<Character> enm = dict.keys(); enm.hasMoreElements();) {
                char nextElement = enm.nextElement();
                //System.out.println(nextElement);
                if (b == nextElement) {
                    int update = dict.get(nextElement);
                    dict.put(b, update+1);
                    break;
                }
            }
        }
        this.composition = dict;
        this.length = sequence.length();
    }

    /**
     * calculates the net charge of a peptide/amino acid sequence
     * @param ph ph value on which the calculation is based
     * @return net charge of the peptide is returned
     */
    public double netCharge(double ph) {
      double C = 8.33;
      double D = 3.86;
      double E = 4.25;
      double H = 6.0;
      double K = 10.53;
      double R = 12.48;
      double Y = 10.07;
      double N_Term =  9.69;
      double C_Term = 2.34;
      Dictionary<Character, Integer> dict = this.composition;
      double left = ((dict.get('R') *(Math.pow(10,R))/(Math.pow(10,ph) + Math.pow(10,R))) +
              (dict.get('K') *(Math.pow(10,K))/(Math.pow(10,ph) + Math.pow(10,K))) +
              (dict.get('H') *(Math.pow(10,H))/(Math.pow(10,ph) + Math.pow(10,H))) +
              (1*(Math.pow(10,N_Term))/(Math.pow(10,ph) + Math.pow(10,N_Term))));
      double right = ((dict.get('D') *(Math.pow(10,ph))/(Math.pow(10,ph) + Math.pow(10,D))) +
              (dict.get('E') *(Math.pow(10,ph))/(Math.pow(10,ph) + Math.pow(10,E))) +
              (dict.get('C') *(Math.pow(10,ph))/(Math.pow(10,ph) + Math.pow(10,C))) +
              (dict.get('Y') *(Math.pow(10,ph))/(Math.pow(10,ph) + Math.pow(10,Y))) +
              (1*(Math.pow(10,ph))/(Math.pow(10,ph) + Math.pow(10,C_Term))));
        double Nettoladung = left - right;
        Nettoladung = Math.round(Nettoladung * 100);
        Nettoladung = Nettoladung / 100;
        this.netCharge = Nettoladung;
        return Nettoladung;
    }

    /**
     * private methode which is only used within the method 'isoelectricPoint'
     * this way (because of the specific structure of above-mentioned function) the net charge of the object is not
     * influenced and should stay correct that way
     * @param ph ph value on which the calculation is based
     * @return net charge of the peptide is returned
     */
    private double netcharge(double ph) {
        double C = 8.33;
        double D = 3.86;
        double E = 4.25;
        double H = 6.0;
        double K = 10.53;
        double R = 12.48;
        double Y = 10.07;
        double N_Term =  9.69;
        double C_Term = 2.34;
        Dictionary<Character, Integer> dict = this.composition;
        double left = ((dict.get('R') *(Math.pow(10,R))/(Math.pow(10,ph) + Math.pow(10,R))) +
                (dict.get('K') *(Math.pow(10,K))/(Math.pow(10,ph) + Math.pow(10,K))) +
                (dict.get('H') *(Math.pow(10,H))/(Math.pow(10,ph) + Math.pow(10,H))) +
                (1*(Math.pow(10,N_Term))/(Math.pow(10,ph) + Math.pow(10,N_Term))));
        double right = ((dict.get('D') *(Math.pow(10,ph))/(Math.pow(10,ph) + Math.pow(10,D))) +
                (dict.get('E') *(Math.pow(10,ph))/(Math.pow(10,ph) + Math.pow(10,E))) +
                (dict.get('C') *(Math.pow(10,ph))/(Math.pow(10,ph) + Math.pow(10,C))) +
                (dict.get('Y') *(Math.pow(10,ph))/(Math.pow(10,ph) + Math.pow(10,Y))) +
                (1*(Math.pow(10,ph))/(Math.pow(10,ph) + Math.pow(10,C_Term))));
        //Nettoladung = Math.round(Nettoladung * 100);
        //Nettoladung = Nettoladung / 100;
        return left - right;
    }

    /**
     * method calculates the iso electric point of a peptide
     * is the pH at which a molecule carries no net electrical charge or is electrically neutral in the statistical mean
     * depending on the net charge, the ph which is used to calculate is either increased (if) or decreased (else) to
     * finally find the ph at which the net charge is zero
     *
     * @return iso electric point is returned
     */
    public double isoelectricPoint() {
        double pHMin = 0;
        double pHMax = 14;
        double pH = (pHMin + pHMax) / 2; // Initial pH guess
        double netCharge = netcharge(pH);

        while (Math.abs(netCharge) >= 0.001) { // user might want to change the tolerance depending on their facilities
            if (netCharge > 0) {
                pHMin = pH;
            } else {
                pHMax = pH;
            }

            pH = (pHMin + pHMax) / 2;
            netCharge = netcharge(pH);
        }
        pH = Math.round(pH * 100);
        pH = pH / 100;
        this.isoelectricPoint = pH;
        return pH;
        //System.out.println(this.Isoelektrischer_Punkt);
    }
    @Override
    public double molecularWeight() {
    return 0;
    }

    @Override
    public double gcContent() {
    return 0;
    }

    @Override
    public void meltingTemp() {
    //empty
    }

    public String nucleotoamino() {return null;}

    /**
     * method to print basic information for a certain fasta object
     *
     * @return basic information about fasta object
     */
    public String getInformation() {
        return ("This Fasta object of type " + this.sequencetype + " with the sequence " + this.sequence + " and the ID " + this.ID +
                " has a length of " + this.length + " and is composed of " + this.composition + " . Further details: net-charge:"
                + this.netCharge + " ,iso electric point: " + this.isoelectricPoint);
    }

}