/*  APBSmem
    Portions Copyright (c) 2010 Michael Grabe
    Portions Copyright (c) 2010-2013 University of Pittsburgh
    Copyright (c) 2014 The Regents of the University of California
    Distributed under The MIT License.  See LICENSE.txt
*/

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package apbsmem;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.Scanner;


/**
 *
 * @author keithc
 */
public class ChargesOff {

    File infile, outfile;
    private Scanner fScan;
    private PrintWriter out;

    /**
     *
     * @param infile
     * @param outfile
     */
    public ChargesOff(String infile, String outfile) {
        this.infile = new File(infile);
        this.outfile = new File(outfile);
    }

    /**
     *
     * @param args argument 1 should be the input PQR file, argument 2 is the output file name
     */
    public static void main(String[] args) {
        ChargesOff co = new ChargesOff(args[0], args[1]);
        co.turnoffcharges();
    }

    /**
     * Read a PQR file, set all atoms to zero charge and save as a new file
     */
    public void turnoffcharges() {

        if(!infile.exists()) {
            System.out.println("ChargesOff error: could not locate pqr file.");
            return;
        }

        try {
            fScan = new Scanner(new FileInputStream(infile));
        } catch (Exception e1) {
            System.out.println(e1.getMessage());
        }

        try {
            out = new PrintWriter(new FileOutputStream(outfile));
        }
        catch (Exception e1) {
               System.out.println("SaveToFile exception: " + e1.toString());
        }

        System.out.println("Turning charges off on input file: " + infile.getName());
        System.out.println("Saving as: " + outfile.getName());

        String nl;

        while (fScan.hasNextLine()) {
            nl = fScan.nextLine();
            if(nl.startsWith("ATOM")) {
                out.println(nl.substring(0, 54) + "  0.0000" + nl.substring(62));
            }
        }

        fScan.close();
        out.close();
    }
}


