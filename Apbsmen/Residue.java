/*  APBSmem
    Portions Copyright (c) 2010 Michael Grabe
    Portions Copyright (c) 2010-2013 University of Pittsburgh
    Copyright (c) 2014 The Regents of the University of California
    Distributed under The MIT License.  See LICENSE.txt
*/

package apbsmem;

import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;

import java.util.List;
import java.util.ArrayList;
import java.util.Arrays;

import java.util.HashMap;

/**
 * Read residue info from a file, extract residues from a file,
 * or modify the charge state of a residue.
 *
 */
public class Residue {
    private File inputFile;
    private String[] residueName;
    private Integer[] residueNumber;
    private Integer[] residueOffset;

    public String[] getResidueName() {
        return residueName;
    }
    public Integer[] getResidueNumber() {
        return residueNumber;
    }
    public Integer[] getResidueOffset() {
        return residueOffset;
    }

    public Residue(String inputFilename) {
        inputFile = new File(inputFilename);

        if (!inputFile.exists()) {
            System.out.println("Residue:  The pqr file " + inputFilename + " does not exist.");
            return;
        }
    }

    /* Read the input pqr file to get a list of the residue names and numbers.
       Residue numbers in the pqr file need not be sequential, so using the order
       in the file is insufficient.

       APBS reports per-atom components numbering atoms sequentially from 0,
       independent of the chain number (field 2) in the pqr file.
     */

    public void readResidueList() throws IOException {
        List<String> resName = new ArrayList<String>();
        List<Integer> resNum = new ArrayList<Integer>();
        List<Integer> resOff = new ArrayList<Integer>();
        int i = 0;

        BufferedReader in = new BufferedReader(new FileReader(inputFile));

        String line;
        String lastName = "";
        int lastNum = -1;

        while ((line = in.readLine()) != null) {

            if (!(line.startsWith("ATOM")|| line.startsWith("HETATM")))
                continue;

            int num = Integer.parseInt(line.substring(21,26).trim());
            if (num != lastNum) {
                resNum.add(lastNum = num);
                resName.add(lastName = line.substring(17,20).trim());
                resOff.add(i);
            }
            i++;
        }

        residueName = resName.toArray(new String[resName.size()]);
        residueNumber = resNum.toArray(new Integer[resNum.size()]);
        residueOffset = resOff.toArray(new Integer[resOff.size()]);

        in.close();
    }

    private String atomName(String line) {
        return line.substring(11,17);
    }

    /* An individual residue in a pqr file consists of a bunch of atoms,
       one per line, in somewhat random order, with each atom being
       identified by a unique name.

       We need to apply the partial charges on the atoms of one residue
       to the atoms of a second residue.  A simple way to do this
       is to create a mapping between atom names and partial charges:

       residue1 = {"C" -> 0.40, "O" -> -0.55, "OD1" -> -0.1, ...}

       Then for each atom in residue 2, look up the partial charge
       in the residue1 map by residue name.

       For our purposes, it make more sense to map the atom name to
       the entire matching line:
       residue1 = {"N" -> "ATOM     38  N   ASP     4      36.757  40.599  28.876 -0.4000 1.5000",
                   ...};

       residueHash(r) returns that map for residue r.
    */

    public HashMap<String,String> residueHash(int r) {
        HashMap<String,String> h = new HashMap<String,String>();
        try {
            BufferedReader in = new BufferedReader(new FileReader(inputFile));
            String line;

            while ((line = in.readLine()) != null) {

                if (!(line.startsWith("ATOM")|| line.startsWith("HETATM")))
                    continue;

                int res = Integer.parseInt(line.substring(21,26).trim());
                if (res != r)
                    continue;

                h.put(atomName(line),line);

                while ((line = in.readLine()) != null) {
                    res = Integer.parseInt(line.substring(21,26).trim());
                    if (res != r)
                        break;
                    h.put(atomName(line),line);
                }
                break;
            }

            in.close();

        } catch (Exception e) {
            System.out.println("Error reading pqr file " + inputFile.getName() + ":  " + e.getMessage());
        }

        return h;
    }

    public void printResidues() {
        for (int i = 0; i < residueName.length; i++) {
            System.out.println("" + i + "  " + residueNumber[i] + "  " + residueName[i]);
        }
    }

    public String writeResidue(int r) {
        String ifn = inputFile.getName();

        String ofn = "";
        int lastDot = ifn.lastIndexOf(".");
        if (lastDot < 0) {  // no "." found in input file name
            ofn = ifn + "." + r;
        } else {
            ofn = ifn.substring(0,lastDot+1) + r + ifn.substring(lastDot);
        }

        return writeResidue(ifn + "." + r, r);
    }

    public String writeResidue(String outputFilename, int r) {
        try {
            BufferedReader in = new BufferedReader(new FileReader(inputFile));
            String line;

            while ((line = in.readLine()) != null) {

                if (!(line.startsWith("ATOM")|| line.startsWith("HETATM")))
                    continue;

                int res = Integer.parseInt(line.substring(21,26).trim());
                if (res != r)
                    continue;

                BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFilename)));
                out.write("REMARK  Residue number " + r + " extracted from file " + inputFile.getName());
                out.newLine();
                out.write(line);
                out.newLine();

                while ((line = in.readLine()) != null) {
                    res = Integer.parseInt(line.substring(21,26).trim());
                    if (res != r)
                        break;
                    out.write(line);
                    out.newLine();
                }

                out.close();
                break;
            }

            in.close();

        } catch (Exception e) {
            System.out.println("Error reading pqr file " + inputFile.getName() + ":  " + e.getMessage());
        }

        return outputFilename;
    }

    public String writeResidueFromModel(int r, HashMap<String,String> model) {
        String ifn = inputFile.getName();

        String ofn = "";
        int lastDot = ifn.lastIndexOf(".");
        if (lastDot < 0) {  // no "." found in input file name
            ofn = ifn + "." + r;
        } else {
            ofn = ifn.substring(0,lastDot+1) + r + ifn.substring(lastDot);
        }

        return writeResidueFromModel(ofn, r, model);
    }

    private String applyModel(String line, HashMap<String,String> model) {
        String name = atomName(line);
        String modline = model.get(name);
        String newline = line.substring(0,54);

        // if the atom doesn't exist in the model, zero the partial charge and radius
        if (null == modline) {
            newline += "  0.0000 0.0000";
        } else {  // otherwise, use the new partial charge and radius, and remove it so we don't reuse it
            newline += modline.substring(54);
            model.remove(name);
        }

        return newline;
    }


    private String writeResidueFromModel(String outputFilename, int r, HashMap<String,String> model) {
        try {
            BufferedReader in = new BufferedReader(new FileReader(inputFile));
            String line;

            while ((line = in.readLine()) != null) {

                if (!(line.startsWith("ATOM")|| line.startsWith("HETATM")))
                    continue;

                int res = Integer.parseInt(line.substring(21,26).trim());
                if (res != r)
                    continue;

                BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFilename)));
                out.write("REMARK  Residue number " + r + " extracted from file " + inputFile.getName() + " with charge model from another pqr file.");
                out.newLine();
                out.write(applyModel(line,model));
                out.newLine();

                while ((line = in.readLine()) != null) {
                    res = Integer.parseInt(line.substring(21,26).trim());
                    if (res != r)
                        break;
                    out.write(applyModel(line,model));
                    out.newLine();
                }

                out.close();
                break;
            }

            in.close();

        } catch (Exception e) {
            System.out.println("Error reading pqr file " + inputFile.getName() + ":  " + e.getMessage());
        }

        return outputFilename;
    }

    private String writeFileFromModel(Residue model, String[] residues) {
        String ifn = inputFile.getName();

        String ofn = "";
        int lastDot = ifn.lastIndexOf(".");
        if (lastDot < 0) {  // no "." found in input file name
            ofn = ifn + ".m";
        } else {
            ofn = ifn.substring(0,lastDot+1) + "m" + ifn.substring(lastDot);
        }

        int[] resses = new int[residues.length];
        for (int i = 0; i < residues.length; i++) {
            resses[i] = Integer.parseInt(residues[i]);
        }
        Arrays.sort(resses);

        return writeFileFromModel(ofn, model, resses);
    }

    private String writeFileFromModel(String outputFilename, Residue model, int[] residues) {
        try {
            BufferedReader in = new BufferedReader(new FileReader(inputFile));
            BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFilename)));
            String line;

            out.write("REMARK  Extracted from file " + inputFile.getName() +
                      " with charge model from " + outputFilename);
            out.newLine();

            int i = 0;
            while ((line = in.readLine()) != null) {

                if (!(line.startsWith("ATOM")|| line.startsWith("HETATM"))) {
                    out.write(line);
                    out.newLine();
                    continue;
                }

                int res = Integer.parseInt(line.substring(21,26).trim());

                if (i == residues.length) {
                    out.write(line);
                    out.newLine();
                    continue;
                }

                if (res > residues[i]) {
                    i++;
                    if (i == residues.length) {
                        out.write(line);
                        out.newLine();
                        continue;
                    }
                }

                if (res != residues[i]) {
                    out.write(line);
                    out.newLine();
                    continue;
                }

                HashMap<String,String> mod = model.residueHash(res);
                out.write(applyModel(line,mod));
                out.newLine();

                while ((line = in.readLine()) != null) {
                    if (!(line.startsWith("ATOM")|| line.startsWith("HETATM"))) {
                        out.write(line);
                        out.newLine();
                        continue;
                    }
                    res = Integer.parseInt(line.substring(21,26).trim());
                    if (res != residues[i]) {
                        i++;
                        out.write(line);
                        out.newLine();
                        break;
                    }
                    out.write(applyModel(line,mod));
                    out.newLine();
                }
            }

            in.close();
            out.close();

        } catch (Exception e) {
            System.out.println("Error reading pqr file " + inputFile.getName() + ":  " + e.getMessage());
        }

        return outputFilename;
    }


    public static void main(String[] args) {
        Residue r = new Residue(args[0]);
        try {
            r.readResidueList();
        } catch (IOException e) {
            System.out.println("Error reading " + args[0]);
            return;
        }

        // See if args[1] points to an existing file
        // If so, apply the charge model from residues in args[1] to the residues in arg[0],
        // zeroing the partial charges of any atoms in residues in args[0] that do not exits in args[1]
        if (args.length > 1) {
            if ((new File(args[1])).exists()) {
                Residue m = new Residue(args[1]);
                try {
                    m.readResidueList();
                } catch (IOException e2) {
                    System.out.println("Error reading " + args[1]);
                    return;
                }

                for (int i = 2; i < args.length; i++) {
                    r.writeResidueFromModel(Integer.parseInt(args[i]),
                                            m.residueHash(Integer.parseInt(args[i])));
                }

                r.writeFileFromModel(m,Arrays.copyOfRange(args,2,args.length));

            } else {
                for (int i = 1; i < args.length; i++) {
                    r.writeResidue(Integer.parseInt(args[i]));
                }
            }
        } else {
            r.printResidues();
        }
    }
}
