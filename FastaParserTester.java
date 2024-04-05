package ProjektOPIS;

import org.yeastrc.proteomics.fasta.FASTAEntry;
import org.yeastrc.proteomics.fasta.FASTAFileParser;
import org.yeastrc.proteomics.fasta.FASTAFileParserFactory;

import java.io.File;
import java.io.IOException;

public class FastaParserTester {
    private void processFASTAFile(String filename) throws IOException {
        FASTAFileParser parser = FASTAFileParserFactory.getInstance().getFASTAFileParser(new File(filename));

        for (FASTAEntry entry = parser.getNextEntry(); entry != null; entry = parser.getNextEntry()) {

            System.out.println("Found " + entry.getHeaders().size() + " headers for this FASTA entry.");
            System.out.println("Found this sequence: " + entry.getSequence());

        }
    }

    public static void main(String[] args) throws Exception {
        FastaParserTester example = new FastaParserTester();
        example.processFASTAFile("/Users/Torben/eclipse-workspace/OPIS/src/Projekt/demo.fasta");
    }
}

