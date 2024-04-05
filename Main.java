package ProjektOPIS;

import org.apache.commons.cli.*;

import java.io.IOException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class Main {
    public static void main(String[] args) throws ParseException, IOException {
        final long timeStart = System.currentTimeMillis();
        // code for reading arguments via command line
        // https://www.youtube.com/watch?v=Za8AQHY-vG8
        CommandLineParser parser = new BasicParser();
        Options options = new Options();

        options.addOption("ip", true, "input pathway for the fasta file");
        options.addOption("op", true, "output pathway for the fasta file");
        options.addOption("sequence_type", true, "sequence type of the fasta");
        options.addOption("threads", true, "define the number of threads");
        options.addOption("h", false, "show help");

        try {
            CommandLine commandLine = parser.parse(options, args);
            System.out.println("Input path: " + commandLine.getOptionValue("ip"));
            System.out.println("Output path: " + commandLine.getOptionValue("op"));
            System.out.println("sequence type: " + commandLine.getOptionValue("sequence_type"));
            System.out.println("You entered the following number of threads: " + commandLine.getOptionValue("threads"));

            if (commandLine.hasOption("h")) {
                HelpFormatter formatter = new HelpFormatter();
                formatter.printHelp("CommandLineParameters", options);

            }
        } catch (ParseException e) {
            throw new RuntimeException(e);
        }
        CommandLine commandLine = parser.parse(options, args); // parsing of cmd line args
        String inputPath = commandLine.getOptionValue("ip"); // setting strings for pathway
        String outputPath = commandLine.getOptionValue("op");
        String sequenceType = commandLine.getOptionValue("sequence_type").toUpperCase(); // use command line parser to get sequence type
        sequencetype typ;
        // setting sequence type, otherwise IllegalArgumentException
        try {
            typ = sequencetype.valueOf(sequenceType);
        } catch (IllegalArgumentException e) {
            throw new IllegalSequenceTypeException("Valid Argument types are DNA, RNA and PEPTIDE");
        }
        FastaAdministration administration = new FastaAdministration();
        int nProzessoren = Runtime.getRuntime().availableProcessors(); // get available processors
        double nThreads;
        // if no argument is given, threads are set to 3/4 of max available
        if (!commandLine.hasOption("threads")) {
            nThreads = nProzessoren * 0.75;

        } else {
            nThreads = Integer.parseInt(commandLine.getOptionValue("threads"));
        }
        // Multi-threading
        if (nProzessoren >= nThreads && nThreads % 2 == 0 || nThreads == 1) {
            ExecutorService executor = Executors.newFixedThreadPool((int) nThreads); // new executor is set up with given number of threads
            String[] fastaFileChunks = administration.splitFastaFile(inputPath, (int) nThreads); // string array which includes the path names of the chunks
            // number of chunks equal to number of threads
            // chunks are given to the multiple executors to be processed
            for (int threadIndex = 0; threadIndex < nThreads; threadIndex++) {
                String fastaChunk = fastaFileChunks[threadIndex];
                executor.execute(() -> {
                    administration.ReadFasta(fastaChunk, typ); // fasta files read and added to the administration
                });
            }
            executor.shutdown(); // shut down executor service
        } else {
            throw new IllegalArgumentException("This device does not provide the number of threads you entered! It must be an even number and less than the max amount available!"); // Exception if unreasonable number of threads is given
        }
        final long timeEnd = System.currentTimeMillis();
        administration.OutputFasta(inputPath, outputPath);
        System.out.println("The program took " + (timeEnd - timeStart) + " milliseconds."); // show how long the program took
    }
}


