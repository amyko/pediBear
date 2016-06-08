package test;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import dataStructures.Pedigree;
import utility.DataParser;

public class TestJays {

	public static void main(String[] args) throws IOException{
	
	String inPath = System.getProperty("user.home") + "/Google Drive/Research/pediBear/data/jays/pedigree.txt";
	String outPath = System.getProperty("user.home") + "/Google Drive/Research/pediBear/data/jays/102jays.true";	
	String indexPath = System.getProperty("user.home") + "/Google Drive/Research/pediBear/data/jays/index2name";	
	int numIndiv = 102;
	
	//get map 
	Map<String, Integer> name2Index = new HashMap<String, Integer>();
	BufferedReader reader = DataParser.openReader(indexPath);
	reader.readLine();
	
	//exclude set
	Set<Integer> exclude = new HashSet<Integer>();
	int[] bad = new int[]{31,32,39,40,41,98,377,382, 22};
	for(int i : bad) exclude.add(i);
	
	String line;
	while((line = reader.readLine())!=null){
		
		String[] fields = line.split("\t");
		
		name2Index.put(fields[1], Integer.parseInt(fields[0]));
		
	}
	
	new Pedigree(inPath, outPath, numIndiv, name2Index, exclude);
	
	
	
	}

	
	
}
