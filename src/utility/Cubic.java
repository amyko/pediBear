package utility;

import java.util.ArrayList;
import java.util.List;

//******************************************************************************
//
// File:    Cubic.java
// Package: benchmarks.detinfer.pj.edu.rit.numeric
// Unit:    Class benchmarks.detinfer.pj.edu.rit.numeric.Cubic
//
// This Java source file is copyright (C) 2008 by Alan Kaminsky. All rights
// reserved. For further information, contact the author, Alan Kaminsky, at
// ark@cs.rit.edu.
//
// This Java source file is part of the Parallel Java Library ("PJ"). PJ is free
// software; you can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation; either
// version 3 of the License, or (at your option) any later version.
//
// PJ is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the GNU General Public License for more details.
//
// A copy of the GNU General Public License is provided in the file gpl.txt. You
// may also obtain a copy of the GNU General Public License on the World Wide
// Web at http://www.gnu.org/licenses/gpl.html.
//
//******************************************************************************



/**
 * Class Cubic solves for the real roots of a cubic equation with real
 * coefficients. The cubic equation is of the form
 * <P>
 * <I>ax</I><SUP>3</SUP> + <I>bx</I><SUP>2</SUP> + <I>cx</I> + <I>d</I> = 0
 * <P>
 * To solve a cubic equation, construct an instance of class Cubic; call the
 * Cubic object's <TT>solve()</TT> method, passing in the coefficients <I>a</I>,
 * <I>b</I>, <I>c</I>, and <I>d</I>; and obtain the roots from the Cubic
 * object's fields. The number of (real) roots, either 1 or 3, is stored in
 * field <TT>nRoots</TT>. If there is one root, it is stored in field
 * <TT>x1</TT>, and fields <TT>x2</TT> and <TT>x3</TT> are set to NaN. If there
 * are three roots, they are stored in fields <TT>x1</TT>, <TT>x2</TT>, and
 * <TT>x3</TT> in descending order.
 * <P>
 * The same Cubic object may be used to solve several cubic equations. Each time
 * the <TT>solve()</TT> method is called, the solution is stored in the Cubic
 * object's fields.
 * <P>
 * The formulas for the roots of a cubic equation come from:
 * <P>
 * E. Weisstein. "Cubic formula." From <I>MathWorld</I>--A Wolfram Web Resource.
 * <A HREF="http://mathworld.wolfram.com/CubicFormula.html" TARGET="_top">http://mathworld.wolfram.com/CubicFormula.html</A>
 *
 * @author  Alan Kaminsky
 * @version 02-Feb-2008
 */
public class Cubic
	{

// Hidden constants.

	private static final double TWO_PI = 2.0 * Math.PI;
	private static final double FOUR_PI = 4.0 * Math.PI;
	
	//for machine precision correction
	private static final double dTol = 1e-6;
	private static final double rootTol = 1e-8;
	
// Exported fields.

	public List<Double> roots;

// Exported constructors.

	/**
	 * Construct a new Cubic object.
	 */
	public Cubic(){}


	public void solve(double a,double b,double c, double d){

		roots = new ArrayList<Double>();
		
		// Verify preconditions.
		if (a == 0.0){
			solve_quadratic(b,c,d);
		}
		if (d==0.0){//ax^3 + bx^2 + cx = x(ax^2+bx+c)
			roots.add(0d);
			solve_quadratic(a,b,c);
		}

		else{
			// Normalize coefficients.
			double denom = a;
			a = b/denom;
			b = c/denom;
			c = d/denom;

			// Commence solution.
			double a_over_3 = a / 3.0;
			double Q = (3*b - a*a) / 9.0;
			double Q_CUBE = Q*Q*Q;
			double R = (9*a*b - 27*c - 2*a*a*a) / 54.0;
			double R_SQR = R*R;
			double D = Q_CUBE + R_SQR;
			
			if (D>0d && D<dTol) D = 0d; //set D=0 (machine precision)

			if (D < 0.0){
				// Three unequal real roots.
				double theta = Math.acos (R / Math.sqrt (-Q_CUBE));
				double SQRT_Q = Math.sqrt (-Q);
				double x1 = 2.0 * SQRT_Q * Math.cos (theta/3.0) - a_over_3;
				double x2 = 2.0 * SQRT_Q * Math.cos ((theta+TWO_PI)/3.0) - a_over_3;
				double x3 = 2.0 * SQRT_Q * Math.cos ((theta+FOUR_PI)/3.0) - a_over_3;
				roots.add(x1);
				roots.add(x2);
				roots.add(x3);
			}
			else if (D > 0.0){
				// One real root.
				double SQRT_D = Math.sqrt (D);
				double S = Math.cbrt (R + SQRT_D);
				double T = Math.cbrt (R - SQRT_D);
				double x1 = (S + T) - a_over_3;
				roots.add(x1);
			}
			else{//debugged this part
				// Three real roots, at least two equal.
				double CBRT_R = Math.cbrt (R);
				double x1 = 2*CBRT_R - a_over_3;
				double x2 = -CBRT_R - a_over_3;
				roots.add(x1);
				roots.add(x2);
			}
		}
		
		roundRoots();
		
		}
	
	public void solve_quadratic (double a, double b, double c){
		
		if (a==0.0){
			solve_linear(b,c);
		}
				
		else{
			double disc = Math.pow(b,2) - 4*a*c;
			
			if (disc==0.0){
				double x1 = -b/(2*a);
				roots.add(x1);
			}
			else if (disc > 0.0){
				double x1 = (-b + Math.sqrt(disc))/(2*a);
				double x2 = (-b - Math.sqrt(disc))/(2*a);
				roots.add(x1);
				roots.add(x2);
			}
			else{
				return;
			}
			
		}
		
	}
	
	public void solve_linear (double a, double b){ //ax + b = 0

		if(a==0.0){
			return;
		}
		else{
			double x1 = -b/a;
			roots.add(x1);
		}
	}
	
	private void roundRoots(){ //round to 0 or 1 for exact solutions
		
		for(int i=0; i<roots.size(); i++){
			if(Math.abs(roots.get(i)-1d)<rootTol){	
				roots.set(i,1d);
			}
			else if (Math.abs(roots.get(i))<rootTol){
				roots.set(i,0d);
			}
		}	
	}
	
	}