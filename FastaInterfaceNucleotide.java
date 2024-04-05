package ProjektOPIS;

public interface FastaInterfaceNucleotide {

	void composition();
	/*
	method to calculate the molecular weight
	 */

	double molecularWeight();
	/*
	method to calculate the percentage of GC
	 */
	double gcContent();

	/*
	method to calculate the melting temperature
	 */
	void meltingTemp();

	/*
	method to translate nucleotide sequence to amino acid sequence
	 */
	String nucleotoamino();

}