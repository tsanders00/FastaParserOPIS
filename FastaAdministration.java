package ProjektOPIS;

import org.yeastrc.proteomics.fasta.FASTAEntry;
import org.yeastrc.proteomics.fasta.FASTAFileParser;
import org.yeastrc.proteomics.fasta.FASTAFileParserFactory;
import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.locks.*;



public class FastaAdministration {

	private int numberFastaEntries = 1;
	//Array for fasta objects
	private final AbstractFastaObjekt[] FastaArray = new AbstractFastaObjekt[9999];
	private final Lock lock = new ReentrantLock(); // lock for sync multi threading

	// empty constructor
	public FastaAdministration() {
	}
	//constructor

	/**
	 * method to read a fasta file and create the corresponding objects
	 * every necessary method is invoked automatically while the fasta file is read
	 * locks are now included to support synchronized multi threading
	 * @param Pfadname absolute path of fasta file
	 * @param typ type of sequence (enum)
	 */
	public void ReadFasta(String Pfadname, sequencetype typ) {
		//https://github.com/yeastrc/java-fasta-utils
		//code for reading a fasta file from a pathway passed earlier and yielding the sequence from it
		FASTAEntry entry;

		try {
			FASTAFileParser fastaparser = FASTAFileParserFactory.getInstance().getFASTAFileParser(new File(Pfadname));
			for (entry = fastaparser.getNextEntry(); entry != null; entry = fastaparser.getNextEntry()) {
				switch (typ) {
					case DNA -> {
						lock.lock();
						try {
							FastaArray[numberFastaEntries] = new DNA();
							FastaArray[numberFastaEntries].addSeq(entry.getSequence().toUpperCase());
							FastaArray[numberFastaEntries].addID(numberFastaEntries);
							FastaArray[numberFastaEntries].setSequenceType(typ);
							FastaArray[numberFastaEntries].composition();
							FastaArray[numberFastaEntries].molecularWeight();
							FastaArray[numberFastaEntries].gcContent();
							FastaArray[numberFastaEntries].meltingTemp();
							//FastaArray[AnzahlFasta].nucleotoamino();
							System.out.println(FastaArray[numberFastaEntries].getInformation());
						} finally {
							lock.unlock();}
					}
					case RNA -> {
						lock.lock();
						try {
							FastaArray[numberFastaEntries] = new RNA();
							FastaArray[numberFastaEntries].addSeq(entry.getSequence().toUpperCase());
							FastaArray[numberFastaEntries].addID(numberFastaEntries);
							FastaArray[numberFastaEntries].setSequenceType(typ);
							FastaArray[numberFastaEntries].composition();
							FastaArray[numberFastaEntries].molecularWeight();
							FastaArray[numberFastaEntries].gcContent();
							FastaArray[numberFastaEntries].meltingTemp();
							System.out.println(FastaArray[numberFastaEntries].getInformation());
						} finally {
							lock.unlock();}
					}
					case PEPTIDE -> {
						lock.lock();
						try {
							FastaArray[numberFastaEntries] = new Peptid();
							FastaArray[numberFastaEntries].addSeq(entry.getSequence().toUpperCase());
							FastaArray[numberFastaEntries].addID(numberFastaEntries);
							FastaArray[numberFastaEntries].setSequenceType(typ);
							FastaArray[numberFastaEntries].composition();
							FastaArray[numberFastaEntries].netCharge(7);
							FastaArray[numberFastaEntries].isoelectricPoint();
							System.out.println(FastaArray[numberFastaEntries].getInformation());
						} finally {
							lock.unlock();}
					}
					case AMBIGUOUS -> throw new AmbiguousStringException("The string type is ambiguous and cannot be processed " +
							"by this programm. Supported types are DNA, RNA and PEPTIDE!");
				}
				lock.lock();
				try {
					numberFastaEntries++;
				} finally {
					lock.unlock();}
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	/**
	 * this method splits the full fasta file into chunks which is necessary for multi threading
	 * the number of chunks is equivalent to the number of threads given by the user
	 * precisely, each fasta entry is written into a file chunk
	 * @param inputFile input fasta file
	 * @param numThreads number of threads given by user
	 * @return is a string array containing the file chunk names
	 * @throws IOException e.g. if no file is found
	 */
	public String[] splitFastaFile(String inputFile, int numThreads) throws IOException {
		List<String> chunkPaths = new ArrayList<>();

		// create output directory
		File outputDir = new File("FastaChunks");
		if (!outputDir.exists()) {
			outputDir.mkdirs();
		}

		// array for buffered writer, each one can write into one file
		// number corresponds to the number of threads
		BufferedWriter[] writers = new BufferedWriter[numThreads];
		for (int i = 0; i < numThreads; i++) {
			String chunkPath = "FastaChunks/chunk" + i + ".fasta"; //chunk number set via index i
			chunkPaths.add(chunkPath); // add path to path list
			writers[i] = new BufferedWriter(new FileWriter(chunkPath));
		}

		try (BufferedReader reader = new BufferedReader(new FileReader(inputFile))) {
			String line;
			int currentThread = 0;
			int lineCount = 0;
			StringBuilder entryBuilder = new StringBuilder();

			while ((line = reader.readLine()) != null) {
				// each entry consists of a header and a sequence which are appended here
				entryBuilder.append(line);
				entryBuilder.append("\n");
				lineCount++;
				// when header and seq are appended, they are written into a chunk
				if (lineCount == 2) {
					// entry ready to be added to the file
					writers[currentThread].write(entryBuilder.toString());
					entryBuilder.setLength(0); // reset string builder
					lineCount = 0; // reset line count

					currentThread = (currentThread + 1) % numThreads;
				}
			}
		}

		// close buffered writer
		for (BufferedWriter writer : writers) {
			if (writer != null) {
				writer.close();
			}
		}

		return chunkPaths.toArray(new String[0]);
	}

	/**
	 * returns the fasta array containing all objects
	 * @return FastaArray
	 */
	public AbstractFastaObjekt[] getFastaArray () {
		return FastaArray;
	}

	/**
	 * method to write all fasta objects into an output file
	 * @param inputPath path for input file to extract file name
	 * @param outputPath path for output file
	 * @throws IOException if no path or file is found
	 */
	public void OutputFasta (String inputPath, String outputPath) throws IOException {
		File inputFile = new File(inputPath); //create file of inputPath to be able to get name of file
		String outputDir;
		// output dir either set through cmd line or where the program is located
		if (outputPath != null) {
			outputDir = outputPath;
		} else {
			outputDir = System.getProperty("user.dir");
		}
		String inputFileName = inputFile.getName(); // set name of input file
		String outputFileName = outputDir + "/" + modifyFileName(inputFileName); // set name of output file
		lock.lock();
		try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFileName))) {
			for (int i = 1; i < numberFastaEntries; i++) {
				AbstractFastaObjekt fastaObject = FastaArray[i];
				String s = fastaObject.getInformation() + "\n"; // create String with information for object i
				writer.write(s); // write String into file
			}
			writer.flush();
		} finally {
			lock.unlock();}
		}

	/**
	 * method to extract the file name from input path
	 * @param inputFileName receives the input fasta file name
	 * @return returns the modified output file name
	 */
	private String modifyFileName(String inputFileName) {
		int dotIndex = inputFileName.lastIndexOf(".");
		String fileNameWithoutExtension = inputFileName.substring(0, dotIndex); // file name
		String fileExtension = inputFileName.substring(dotIndex); // file format
		return fileNameWithoutExtension + "_output" + fileExtension;
	}
	/**
	 * returns the number of fasta entries
	 * @return AnzahlFasta
	 */
	public int getNumberOfEntries() {
		return this.numberFastaEntries;
	}

}
