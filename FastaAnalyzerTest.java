package ProjektOPIS;
import org.junit.Test;
import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 *
 */
public class FastaAnalyzerTest {

    @Test
    public void testDNA() {
        FastaAdministration vw = new FastaAdministration();
        sequencetype typ = sequencetype.valueOf("DNA");
        vw.ReadFasta("/Users/Torben/eclipse-workspace/OPIS/src/ProjektOPIS/testDNA.fasta", typ);
        AbstractFastaObjekt[] arr = vw.getFastaArray();
        //https://molbiotools.com/dnacalculator.php
        //values are rounded
        assertEquals(48, arr[0].gcContent());
        assertEquals(15452, arr[0].molecularWeight());

    }

    @Test
    public void testRNA() {
        FastaAdministration vw = new FastaAdministration();
        sequencetype typ = sequencetype.valueOf("RNA");
        vw.ReadFasta("/Users/Torben/eclipse-workspace/OPIS/src/ProjektOPIS/testRNA.fasta", typ);
        AbstractFastaObjekt[] arr = vw.getFastaArray();
        //https://molbiotools.com/dnacalculator.php
        //values are rounded
        assertEquals(49, arr[0].gcContent());
        assertEquals(22610, arr[0].molecularWeight());
        //usual error is arround 1 to 5 %
        //acceptable?

    }

    @Test
    public void testPeptid() {
        FastaAdministration vw = new FastaAdministration();
        sequencetype typ = sequencetype.valueOf("PEPTIDE");
        vw.ReadFasta("/Users/Torben/eclipse-workspace/OPIS/src/ProjektOPIS/testPeptid.fasta", typ);
        AbstractFastaObjekt[] arr = vw.getFastaArray();
        //https://www.bachem.com/knowledge-center/peptide-calculator/
        assertEquals(7.24,arr[0].netCharge(7));
        assertEquals(10.1 , arr[0].isoelectricPoint());
        //result should be 10.10 but is calculated as 10.09
        //calculation still delivers good and precise results
        //rounding error not really relevant
    }

    @Test
    public void testErrorMessage() {
        try {
            // Call the method that is expected to throw an exception
            FastaAdministration vw = new FastaAdministration();
            sequencetype typ = sequencetype.valueOf("AMBIGUOUS");
            vw.ReadFasta("/Users/Torben/eclipse-workspace/OPIS/src/ProjektOPIS/testDNA.fasta", typ);
        } catch (Exception AmbiguousStringException) {
            // Check if the error message matches the expected value
            assertEquals("The string type is ambiguous and cannot be processed by this programm. Supported types are DNA, RNA and PEPTIDES!", AmbiguousStringException.getMessage());
        }
    }
}
