package utility;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.FileReader;
import java.io.BufferedReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

public class DataParser {	

	public static BufferedReader openReader(String path) throws FileNotFoundException{
		// open file reader
		FileReader fr = new FileReader(path);
		BufferedReader textReader = new BufferedReader(fr);
		return textReader;
	}
	
	public static PrintWriter openWriter(String path) throws IOException{
		//open writers
		FileWriter write = new FileWriter(path);
		PrintWriter print_ln = new PrintWriter(write);
		return print_ln;
	}
	
	public static int countLines(String path) throws IOException{
		BufferedReader reader = openReader(path);
		int count = 0;
		while (reader.readLine() != null) count++;
		reader.close();
		return count;
	}
	
	//concatenate files given by the list of paths at specified columns
	public static void concatFiles(String[] paths, String outPath, int[][] cols) throws IOException{
		//assumes that the line nums are equal //fix this later
		
		//open readers & writer
		int numFiles = paths.length;
		List<BufferedReader> readers = new ArrayList<BufferedReader>(numFiles);
		for(String p : paths){
			readers.add(openReader(p));
		}
		PrintWriter writer = openWriter(outPath);
		
		
		// write to file
		boolean endOfLine = false;
		while(true){
			//read lines
			String[] lines = new String[numFiles];
			for (int i=0; i<numFiles; i++){
				lines[i] = readers.get(i).readLine();
				if (lines[i]==null) endOfLine=true;
			}
		
			//check end of line
			if(endOfLine==true) break;
			
			
			//write to file
			String toWrite = "";
			for (int i=0; i<numFiles; i++){
				String[] fields = lines[i].split("\\s");
				
				//write all columns
				if(cols[i].length==0){
					for (int j=0; j<fields.length; j++){
						toWrite += fields[j]+"\t";
					}
				}
				
				//write select columns
				else{
					for (int j=0; j<cols[i].length; j++){
						toWrite += fields[cols[i][j]]+"\t";
					}
				}
			}
			writer.write(toWrite + "\n");	
		}
		
		
		//close files
		for (BufferedReader r : readers){
			r.close();
		}
		writer.close();
		
	}

	public static void concatFilesSpace(String[] paths, String outPath, int[][] cols) throws IOException{
		//assumes that the line nums are equal //fix this later
		
		//open readers & writer
		int numFiles = paths.length;
		List<BufferedReader> readers = new ArrayList<BufferedReader>(numFiles);
		for(String p : paths){
			readers.add(openReader(p));
		}
		PrintWriter writer = openWriter(outPath);
		
		
		// write to file
		boolean endOfLine = false;
		while(true){
			//read lines
			String[] lines = new String[numFiles];
			for (int i=0; i<numFiles; i++){
				lines[i] = readers.get(i).readLine();
				if (lines[i]==null) endOfLine=true;
			}
		
			//check end of line
			if(endOfLine==true) break;
			
			
			//write to file
			String toWrite = "";
			for (int i=0; i<numFiles; i++){
				String[] fields = lines[i].split("\\s");
				
				//write all columns
				if(cols[i].length==0){
					for (int j=0; j<fields.length; j++){
						toWrite += fields[j]+" ";
					}
				}
				
				//write select columns
				else{
					for (int j=0; j<cols[i].length; j++){
						toWrite += fields[cols[i][j]]+" ";
					}
				}
			}
			writer.write(toWrite + "\n");	
		}
		
		
		//close files
		for (BufferedReader r : readers){
			r.close();
		}
		writer.close();
		
	}
	
	public static void writeMatrixToFile(PrintWriter writer, double[][] data) throws IOException{
		
		//write data
		for (int i=0; i<data.length; i++){
			for(int j=0; j<data[0].length; j++){
				writer.write(String.format("%f\t",data[i][j]));
			}
			writer.write("\n");
		}
		
	}
	
	public static void writeStringArray(PrintWriter writer, String[] toWrite){
		for(String item : toWrite){
			writer.write(item+"\t");
		}
		writer.write("\n");
	}
}