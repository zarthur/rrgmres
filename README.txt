  ***************************************************************************
  * All the software  contained in this library  is protected by copyright. *
  * Permission  to use, copy, modify, and  distribute this software for any *
  * purpose without fee is hereby granted, provided that this entire notice *
  * is included  in all copies  of any software which is or includes a copy *
  * or modification  of this software  and in all copies  of the supporting *
  * documentation for such software.                                        *
  ***************************************************************************
  * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED *
  * WARRANTY. IN NO EVENT, NEITHER  THE AUTHORS, NOR THE PUBLISHER, NOR ANY *
  * MEMBER  OF THE EDITORIAL BOARD OF  THE JOURNAL  "NUMERICAL ALGORITHMS", *
  * NOR ITS EDITOR-IN-CHIEF, BE  LIABLE FOR ANY ERROR  IN THE SOFTWARE, ANY *
  * MISUSE  OF IT  OR ANY DAMAGE ARISING OUT OF ITS USE. THE ENTIRE RISK OF *
  * USING THE SOFTWARE LIES WITH THE PARTY DOING SO.                        *
  ***************************************************************************
  * ANY USE  OF THE SOFTWARE  CONSTITUTES  ACCEPTANCE  OF THE TERMS  OF THE *
  * ABOVE STATEMENT.                                                        *
  ***************************************************************************

   AUTHOR:

    A. Neuman
        DEPARTMENT OF MATHEMATICAL SCIENCES
        KENT STATE UNIVERSITY

    L. Reichel
        DEPARTMENT OF MATHEMATICAL SCIENCES
        KENT STATE UNIVERSITY

    H. Saddok
        LABORATOIRE DE MATHEMATIQUES PURES ET APPLIQUEES
        UNIVERSITE DU LITTORAL
    
   REFERENCE:



   SOFTWARE REVISION:

       Ver 1.0  MAY 2011

   SOFTWARE LANGUAGE:

       MATLAB 7.11

**************************************************************************

Range Restricted GMRES Demo.
 Version 1.0  May 2011.
 Copyright (c) 2011

The installation of Regularization Tools is very simple:

=======
Windows
=======

Step 1. Download the na33 package from Netlib

Step 2. Download the MATLAB package Regularization Tools by Hansen
	from http://www2.imm.dtu.dk/~pch/Regutools/index.html

Step 3. Unzip the na33 package and Regularization Tools to a common
	directory (Use, for example, Winzip to extract all files 
        and store them in RRGDEMO.  The files from Regularization
	Tools will be extracted to REGU)

Step 4. Start MATLAB

Step 5. Add a path to the directories RRGDEMO and REGU

Step 6. Remove the package files


=======
Linux
=======

Step 1. Download the na33 package from Netlib

Step 2. Download the MATLAB package Regularization Tools by Hansen
	from http://www2.imm.dtu.dk/~pch/Regutools/index.html

Step 2. Unzip the packages by using the command
        unzip package_name, where package_name is the file to be
	unzipped.  The directory RRGDEMO will be created by extracting
	the na33 package.  The files from the Regularization Tools
	will be extracted to REGU.) 

Step 3. Start MATLAB

Step 4. Add a path to the directories RRGDEMO and REGU

Step 5. Remove the files with the command rm file_name, where
	file_name is the name of the file to be removed.


The file primer.pdf contains the manual in PDF form.
