package ProjektOPIS;

public abstract class AbstractFastaObjekt implements FastaInterfaceNucleotide,FastaInterfacePeptide {

	/*
	method to return the above information for a fasta object
	 */
	public String getInformation() {
		return null;
	}
	/*
	method to add sequence to an object
	 */
	public void addSeq(String sequence) {
	}
	/*
	method to set an ID
	 */
	public void addID(int ID) {}

	/*
    method to set the sequence_type
     */
	public void setSequenceType(sequencetype typ) {}
}