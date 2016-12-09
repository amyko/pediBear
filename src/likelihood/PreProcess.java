package likelihood;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Map;

import executables.Run;
import utility.DataParser;

// 1) has tped, tfam files, 2) valid allele characters, 3) increasing distance

public class PreProcess {
	
	private static int pA = 7;
	private static int pC = 8;
	private static int pG = 9;
	private static int pT = 10;
	
	
	//check if options are valid and set values
	public static void processOptionfile(String args[]) throws IOException{
		
		//read option file from command line
		if(args.length == 0){
			System.out.println("Missing option file! Exiting program.");
			System.exit(1);
		}
		else if(args.length > 1){
			System.out.println(String.format("Detected more than one command line argument. Reading %s as option file and ignoring the rest.", args[0]));
			System.exit(1);
		}
		
		
		//check that option file can be read
		checkFile(args[0]);
		
		
		
		//check that file name is given
		BufferedReader infile = DataParser.openReader(args[0]);
		String line;
		
		while((line=infile.readLine())!=null){
			
			String[] fields = line.split("\\s#");
			String key = fields[1].split("\\s+")[0].toLowerCase();
			
			//error checking
			if(key.equals("filename") || key.equals("refpopfilename")){
				
				if(fields.length<2){
					System.out.println("No files given. Exiting program.");
					System.exit(1);
				}
				
				
				else if(!checkFile(fields[0]+".tped") || !checkFile(fields[0]+".tfam")){
					System.out.println(String.format("%s does not exist or cannot be opened! Existing program.", fields[0]));
					System.exit(1);
				}
				
				else{		
					if(key.equals("filename")) Run.fileName = fields[0];
					else if(key.equals("refpopfilename")) Run.refPopFileName = fields[0];
					
					//set numIndiv
					Run.numIndiv = DataParser.countLines(Run.fileName+".tfam");
					
				}
				
			}
			
			else if(key.equals("agefilename")){
				
				if(fields.length<2){
					System.out.println("No files given. Exiting program.");
					System.exit(1);
				}
				else if(!checkFile(fields[0])){
					System.out.println(String.format("%s does not exist or cannot be opened! Existing program.", fields[0]));
					System.exit(1);
				}
				
				else{		
					Run.ageFileName = fields[0];
				}
				
				
			}
			
			else if(key.equals("maf")){
			
				if(fields.length < 2){
					System.out.println("MAF not given. Using default value of 0.01");
				}
				else if(!checkNumeric(fields[0], 0, 1)){
					System.out.println("Invalid MAF given. Using default value of 0.01");
				}
				else{
					Run.maf = Double.parseDouble(fields[0]);
				}
			}
			
			
			else if(key.equals("errorrate")){
				
				if(fields.length < 2){
					System.out.println("ErrorRate not given. Using default value of 0.01");
				}
				else if(!checkNumeric(fields[0], 0, 1)){
					System.out.println("Invalid ErrorRate given. Using default value of 0.01");
				}
				else{
					Run.errorRate = Double.parseDouble(fields[0]);
				}
				
			}
			
			else if(key.equals("maxgen")){
				
				if(fields.length < 2){
					System.out.println("ErrorRate not given. Using default value of 5.");
				}
				else if(!checkNumeric(fields[0], 0, 6)){
					System.out.println("Invalid ErrorRate given. Using default value 5.");
				}
				else{
					Run.maxDepth = Integer.parseInt(fields[0]) - 1;
				}
				
			}
			
			else if(key.equals("maxsampledepth")){
				
				if(fields.length < 2){
					System.out.println("sampleGen not given. Using default value of 5");
				}
				else if(!checkNumeric(fields[0], 0, 6)){
					System.out.println("Invalid sampleGen given. Using default value 5");
				}
				else{
					Run.sampleDepth = Integer.parseInt(fields[0]) - 1;
				}
				
			}
			
			
			else if(key.equals("back")){
				
				if(fields.length < 2){
					System.out.println("back not given. Using default value of 0.04");
				}
				else if(!checkNumeric(fields[0], 0, 2)){
					System.out.println("Invalid back given. Using default value 0.04");
				}
				else{
					Run.back = Double.parseDouble(fields[0]);
				}
				
			}
			
			
			else if(key.equals("starttemp")){
				
				if(fields.length < 2){
					System.out.println("startTemp not given. Using default value of 100");
				}
				else if(!checkNumeric(fields[0], 0, 1000)){
					System.out.println("Invalid startTemp given. Using default value 100");
				}
				else{
					Run.startTemp = Double.parseDouble(fields[0]);
				}
				
			}
			
			else if(key.equals("iterpertemp")){
				
				if(fields.length < 2){
					System.out.println("IterPerTemp not given. Using default value of 20000");
				}
				else if(!checkNumeric(fields[0], 0, Double.POSITIVE_INFINITY)){
					System.out.println("Invalid IterPerTemp given. Using default value 20000");
				}
				else{
					Run.iterPerTemp = Integer.parseInt(fields[0]);
				}
				
			}
			
			else if(key.equals("tempfact")){
				
				if(fields.length < 2){
					System.out.println("TempFact not given. Using default value of 1.01");
				}
				else if(!checkNumeric(fields[0], 1, Double.POSITIVE_INFINITY)){
					System.out.println("Invalid TempFact given. Using default value 1.01");
				}
				else{
					Run.tempFact = Double.parseDouble(fields[0]);
				}
				
			}
			
			else if(key.equals("maxiter")){
				
				if(fields.length < 2){
					System.out.println("MaxIter not given. Using default value of 4000000");
				}
				else if(!checkNumeric(fields[0], 1, Double.POSITIVE_INFINITY)){
					System.out.println("Invalid maxIter given. Using default value 4000000");
				}
				else{
					Run.maxIter = Integer.parseInt(fields[0]);
				}
				
			}
			
			
			else if(key.equals("conv")){
				
				if(fields.length < 2){
					System.out.println("conv not given. Using default value of 1");
				}
				else if(!checkNumeric(fields[0], 0, Double.POSITIVE_INFINITY)){
					System.out.println("Invalid conv given. Using default value 1");
				}
				else{
					Run.conv = Double.parseDouble(fields[0]);
				}
				
			}
			
			else if(key.equals("conditiononld")){
				
				if(fields.length < 2){
					System.out.println("conditionOnLD not given. Using default value of 1");
				}
				else if(!fields[0].equals("0") && !fields[0].equals("1")){
					System.out.println("Invalid conditionOnLD given. Using default value 1");
				}
				else{
					Run.conditional = Integer.parseInt(fields[0])==1 ? true : false;
				}
				
			}
			
			else if(key.equals("numrun")){
				
				if(fields.length < 2){
					System.out.println("numRun not given. Using default value of 3");
				}
				else if(!checkNumeric(fields[0],0,100)){
					System.out.println("Invalid numRun given. Using default value 3");
				}
				else{
					Run.numRun = Integer.parseInt(fields[0]);
				}
				
			}
			
			
			else if(key.equals("poissonmean")){
				
				if(fields.length < 2){
					System.out.println("poissonMean not given. Using default value of numIndiv");
				}
				else if(!checkNumeric(fields[0], 0, Double.POSITIVE_INFINITY)){
					System.out.println("Invalid numRun given. Using default value numIndiv");
				}
				else{
					Run.poissonMean = Double.parseDouble(fields[0]);
				}
				
			}
			
			else if(key.equals("numthreads")){
				
				if(fields.length < 2){
					System.out.println("numThreads not given. Using default value of 1");
				}
				else if(!checkNumeric(fields[0], 0, Double.POSITIVE_INFINITY)){
					System.out.println("Invalid numThread given. Using default value 1");
				}
				else{
					Run.numThreads = Integer.parseInt(fields[0]);
				}
				
			}
			
			//option not recognized
			else{
				System.out.println(String.format("Unrecognized option: %s. Exiting program.", key));
				System.exit(1);
			}			
			
		}
		
		infile.close();
		
		
		//check files are given
		if(Run.fileName.equals("")){
			System.out.println("No file given. Exiting program");
			System.exit(1);
		}
		//check refpop file
		if(Run.refPopFileName.equals("")){
			System.out.println("RefPop not given. Setting refPopFileName = fileName");
			Run.refPopFileName = Run.fileName;
		}
		//check depth
		if(Run.sampleDepth > Run.maxDepth){
			System.out.println("Number of generations spanned by samples exceeds max generation. Setting sampleGen = maxGen.");
			Run.sampleDepth = Run.maxDepth;
		}
		//check poisson mean
		if(Run.poissonMean == 0){
			Run.poissonMean = Run.numIndiv;
		}
		
		
		
		
	}
	

	//check that string is a number and that it falls between min/max
	public static boolean checkNumeric(String myString, double min, double max){
		
		try{
			double val = Double.parseDouble(myString);
			if(val > min && val < max) return true;
			return false;
			
		}
		catch(NumberFormatException nfe){
			return false;
		}
		
	}
	
	
	
	//check that file exists & is readable
	public static boolean checkFile(String inPath){
		
		File f = new File(inPath);
		
		if(!f.exists() || !f.isFile() || !f.canRead()){
			return false;
		}
		
		return true;
				
	}
	
	
	public static void checkInputFiles(String testPath, String popPath) throws IOException{
		
		boolean sameFile = testPath.equals(popPath);
		
		//check that both files have the same set of markers
		if(!sameFile)
			sameMarkers(testPath, popPath);
		
		//check tped/tfam files
		checkInput(testPath);
		if(!sameFile)
			checkInput(popPath);
	}
	
	
	//check that test.tped and pop.tped have the same set of SNPs
	public static void sameMarkers(String testPath, String popPath) throws IOException{
		
		BufferedReader infile1 = DataParser.openReader(testPath+".tped");
		BufferedReader infile2 = DataParser.openReader(popPath+".tped");
		String line1;
		String line2;

		
		while((line1=infile1.readLine())!=null && (line2=infile2.readLine())!=null){
		
			String[] fields1 = line1.split("\\s");
			String[] fields2 = line2.split("\\s");
			
			for(int i=0; i<4; i++){
				if(!fields1[i].equals(fields2[i])){
					
					System.out.println(String.format("Sample individuals and reference population do not have the same set of markers. Check that %s and %s have the same marker positions/distances.", testPath, popPath));
					System.exit(1);
					
				}
			}
			
			
		}
		
		infile1.close();
		infile2.close();
		
	}
	
	
	
	
	//check that tped & tfam files are kosher
	public static void checkInput(String testPath) throws IOException{

		
		//get number of individuals
		boolean correctNumIndiv = false;

		
		
		//check line-by-line
		BufferedReader infile = DataParser.openReader(testPath+".tped");
		String line;
		String prevChrom = "NA";
		String currChrom;
		double prevPos = -1;
		double currPos = 0;
		int[] count = new int[4];
		
		while((line=infile.readLine())!=null){
			
			String[] fields = line.split("\\s");
			currChrom = fields[0];
			currPos = Double.parseDouble(fields[2]);
			for(int k=0; k<count.length; k++) count[k] = 0;
			
			//if new chromosome
			if(!currChrom.equals(prevChrom)){
				prevPos = -1;
				
				if(currChrom=="X" || currChrom=="Y"){
					System.out.println(String.format("Sex chromosome detected in %s. Please remove them.", testPath+".tped"));
				}
			}
			
			//check numIndiv matches num genotypes
			if(correctNumIndiv==false){	
				
				int numIndiv = DataParser.countLines(testPath+".tfam");
				
				if((fields.length - 4)/2 != numIndiv){
					System.out.println(String.format("Number of individuals in %s does not match number of genotypes in %s.", testPath+".tfam", testPath+".tped"));	
					System.exit(1);
				}
				else correctNumIndiv = true;
			}
			
			//check distance
			if(currPos==prevPos){	
				System.out.println(String.format("Thre are more than one markers in the same position in %s! Every marker must have a unique position.", testPath+".tped"));
				System.exit(1);
			}
			if(currPos < prevPos) {
				System.out.println(String.format("Detected decreasing distance in %s. Exiting program", testPath+".tped"));
				System.exit(1);
			}
			
			//check biallelic
			for(int i=4; i<fields.length; i++){
				
				String g = fields[i];
				
				if(g.equals("A")) count[0]++;
				else if(g.equals("C")) count[1]++;
				else if(g.equals("G")) count[2]++;
				else if(g.equals("T")) count[3]++;
				//else if(g.equals("0")){}
				else{
					System.out.println(String.format("Illegal allele present in %s. All entires must be either A,C,G, or T. Exiting program.", testPath+".tped"));
					System.exit(1);
				}
				
			}
			
			int nAlleles = 0;
			for(int i=0; i<count.length; i++){
				if(count[i]>0) nAlleles++;
			}
			
			if(nAlleles > 2){
				System.out.println(String.format("Encountered non-biallelic sites in %s. Please filter out non-biallelic sites. Exiting program.", testPath));
				System.exit(1);
			}
			
			
			//update prev
			prevChrom = currChrom;
			prevPos = currPos;
			
			
		}
		
		
	}
	

	
	
	

}
