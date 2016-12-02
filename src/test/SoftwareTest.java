package test;

import java.io.IOException;

import likelihood.LDStreamPedMissing;

public class SoftwareTest {
	
	
	public static void main(String[] args) throws IOException{
		
		String dir = "/Users/kokocakes/Google Drive/Research/pediBear/data/simulations/softwareTest/";
		String fileName = dir + "test";
		String outPath = dir + "test.info";
		int back = 100;
		boolean conditional = true;


		LDStreamPedMissing.writeLdOutfile(fileName, outPath, back, conditional);
	}
	
	
}
