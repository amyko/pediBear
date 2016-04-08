package utility;

//import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import dataStructures.Node;
import dataStructures.Pedigree;
import Unused.PedigreeNode;

public class ArrayUtility {
	
	public static final double EPSILON = 1e-14;
	
	public static double sum(double[] array){
		double toReturn = 0d;
		for( double val : array){
			toReturn += val;
		}
		return toReturn;
	}
	
	public static int drawFrom(double[] dist, Random rgen){
		assert Math.abs(sum(dist) - 1) < EPSILON;
		int idx = -1;
		double p = rgen.nextDouble();
		while(p > 0){
			idx += 1;
			p -= dist[idx];
		}
		return idx;
	}
	
	
	public static void clear(int[] array){
		
		for(int i=0; i< array.length; i++){
			array[i] = 0;
		}
		
	}
	

	
	
	public static List<PedigreeNode> shallowCopyList (List<PedigreeNode> toBeCloned){
		
		List<PedigreeNode> clone = new ArrayList<PedigreeNode>();
		
		for (PedigreeNode item : toBeCloned){
			clone.add(item);
		}
		
		return clone;
	}
	
	
	
	public static int[] getNRandomIdx(int N, int sampleSize, Random rGen){//works

		int[] toReturn = new int[sampleSize];

		int nSampled = 0;
		int i = 0;

		while(nSampled < sampleSize){
			
			if(rGen.nextDouble()*(N-i) < (sampleSize - nSampled)){
				toReturn[nSampled++] = i;
			}
			
			i++;
			
		}
		
		return toReturn;
	}

	
	
	/*
	public static List<PedigreeNode> deepCloneList (List<PedigreeNode> toBeCloned){
		
		List<PedigreeNode> clone = new ArrayList<PedigreeNode>();
		
		for (PedigreeNode item : toBeCloned){
			clone.add((PedigreeNode)deepClone(item));
		}
		
		return clone;
		
	}
	
	public static Object deepClone(Object object) {
		try {
			ByteArrayOutputStream baos = new ByteArrayOutputStream();
		    ObjectOutputStream oos = new ObjectOutputStream(baos);
		    oos.writeObject(object);
		    ByteArrayInputStream bais = new ByteArrayInputStream(baos.toByteArray());
		    ObjectInputStream ois = new ObjectInputStream(bais);
		    return ois.readObject();
		}
		catch (Exception e) {
			e.printStackTrace();
		    return null;
		}
	}
	*/
	
}
