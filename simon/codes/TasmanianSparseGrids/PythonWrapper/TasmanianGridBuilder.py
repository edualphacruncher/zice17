#!/usr/bin/python

import os
import re
import string
import sys

bDebugInfo = False

lTasGridBuilderKnownPoints = [ 
        "clenshaw-curtis",  "clenshaw-curtis-zero",
        "chebychev"  "chebyshev-odd",
        "gauss-legendre", "gauss-legendre-odd",
        "fejer-2", "leja", "leja-odd",
        "max-lebesgue",  "max-lebesgue-odd",  "min-lebesgue", "min-lebesgue-odd",
        "custom-tabulated", "gauss-chebyshev-1",  "gauss-chebyshev-1-odd",
        "gauss-chebyshev-2", "gauss-chebyshev-2-odd",
        "gauss-gegenbauer",  "gauss-gegenbauer-odd",
        "gauss-jacobi", "gauss-laguerre", "gauss-hermite", "gauss-hermite-odd" ]

lTasGridBuilderGlobalRefine = [ 
        "leja", "leja-odd", "max-lebesgue",  "max-lebesgue-odd",  "min-lebesgue", "min-lebesgue-odd" ]

lTasGridBuilderRequireAlpha = [
        "gauss-gegenbauer", "gauss-gegenbauer-odd",
        "gauss-jacobi", "gauss-laguerre", "gauss-hermite", "gauss-hermite-odd" ]

lTasGridBuilderRequireBeta = [
        "gauss-hermite", "gauss-hermite-odd" ]

lTasGridBuilderKnownTypes = [ "polynomial", "polynomial-zero-boundary", "wavelet" ]

lTasGridBuilderKnownRefinementTypes = [ "classic", "family", "direction", "fds" ]

lTasGridBuilderKnownCLICommands = [ "-s", "--silent", "-v", "--verbose", "--reset", "--recover" ]


class TasGridWrapper:
        def __init__( self ):
                pass
                
        def writeMatrix( self, llfMatrix, sFileName ):
                sLines = []
                iRows = len( llfMatrix )
                iCols = len( llfMatrix[0] )
                sS = "{0:1d} {1:1d}\n".format( iRows, iCols )
                sLines.append( sS )
                
                for lRow in llfMatrix:
                        sS = ""
                        for fVal in lRow:
                                sS += "{0:2.17e} ".format( fVal )
                        sS += "\n"
                        sLines.append( sS )
                        
                F = open( sFileName, 'w' )
                for sS in sLines:
                        F.write( sS )
                F.close()
                
        def readMatrix( self, sFileName ):
                F = open( sFileName, 'r' )
                lsLines = F.readlines()
                F.close()
                
                lL = re.split( " ", lsLines[0].strip() )
                
                iRows = int(lL[0])
                iCols = int(lL[1])
                
                llfMatrix = []
                
                for iI in range( len(lsLines) -1 ):
                        sLine = lsLines[iI+1]
                        
                        sLine = string.replace( sLine, "\t", " " )
                        while ( "  " in sLine ):
                                sLine = string.replace( sLine, "  ", " " )
                        sLine = sLine.strip()
                        
                        lL = re.split( " ", sLine )
                        
                        lRow = []
                        for sL in lL:
                                lRow.append( float( sL ) )
                        
                        llfMatrix.append( lRow )
                        
                return llfMatrix
                
        def writeInputFile( self, lfVector, sFileName ):
                sLines = []
                for fVal in lfVector:
                        sLines.append( "{0:2.17e}\n".format( fVal ) )
                F = open( sFileName, 'w' )
                for sS in sLines:
                        F.write( sS )
                F.close()
                
        def readOutputFile( self, sFileName ):
                F = open( sFileName, 'r' )
                lsLines = F.readlines()
                F.close()
                
                lfVector = []
                for sLine in lsLines:
                        lfVector.append( float( sLine ) )
                        
                return lfVector
                        
                
        def readDescription( self, bSilent ):
                # reads ProblemDescription.in, returns true if successful
                if ( not os.path.exists('ProblemDescription.in') ):
                        print(" ERROR: Cannot find ProblemDescription.in")
                        return False
                        
                # read the file
                F = open( 'ProblemDescription.in', 'r' )
                lPDLines = F.readlines()
                F.close()
                
                # parse the input
                iCountLines = 1
                for sRawLine in lPDLines:
                        #process each line
                        lComments = re.split( '#', sRawLine ) # strip comments
                        sRawLine = lComments[0]
                        
                        sLine = sRawLine.strip("\n\r\t")
                        sLine = sLine.strip()
                        
                        if ( len( sLine ) > 0 ):
                                # if the line is not empty or just comment, process the command
                                lCommand = re.split( ':', sLine )
                                
                                if ( len( lCommand ) != 2 ):
                                        print(" ERROR: wroing format, commands should be one per line and name and value should be separated by colon")
                                        print("        at line {0:3d}".format( iCountLines ))
                                        return False
                                
                                lCommand[0] = lCommand[0].strip()
                                lCommand[1] = lCommand[1].strip()
                                
                                if ( (len(lCommand[0]) == 0) or (len(lCommand[1]) == 0) ):
                                        print(" ERROR: missing command or value")
                                        print("        at line {0:3d}".format( iCountLines ))
                                        return False
                                        
                                lCommand[0] = string.replace( lCommand[0], "\t", " " )
                                while ( "  " in lCommand[0] ):
                                        lCommand[0] = string.replace( lCommand[0], "  ", " " )
                                lCommand[0] = lCommand[0].strip()
                                
                                if ( "model executable file" in lCommand[0].lower() ):
                                        self.sModelExecutableFile = lCommand[1]
                                elif ( "model input file" in lCommand[0].lower() ):
                                        self.sModelInputFile = lCommand[1]
                                elif ( "model output file" in lCommand[0].lower() ):
                                        self.sModelOutputFile = lCommand[1]
                                elif ( "tasgrid executable" in lCommand[0].lower() ):
                                        self.sTasgridExec = lCommand[1]
                                elif ( "project name" in lCommand[0].lower() ):
                                        self.sProjectName = lCommand[1]
                                elif ( "number of inputs" in lCommand[0].lower() ):
                                        self.iNumDimensions = int( lCommand[1] )
                                        if ( (self.iNumDimensions < 1) ):
                                                print(" ERROR: inputs should be a positive integer")
                                                print("        at line {0:3d}".format( iCountLines ))
                                                return False
                                elif ( "number of outputs" in lCommand[0].lower() ):
                                        self.iNumOutputs = int( lCommand[1] )
                                        if ( (self.iNumOutputs < 1) ):
                                                print(" ERROR: inputs should be a positive integer")
                                                print("        at line {0:3d}".format( iCountLines ))
                                                return False
                                elif ( "domain transform a" in lCommand[0].lower() ):
                                        lCommand[1] = string.replace( lCommand[1], "\t", " " )
                                        while ( "  " in lCommand[1] ):
                                                lCommand[1] = string.replace( lCommand[1], "  ", " " )
                                        lVals = re.split( " ", lCommand[1] )
                                        self.lfTransformA = []
                                        for sVal in lVals:
                                                self.lfTransformA.append( float( sVal ) )
                                        #print(self.lfTransformA)
                                elif ( "domain transform b" in lCommand[0].lower() ):
                                        lCommand[1] = string.replace( lCommand[1], "\t", " " )
                                        while ( "  " in lCommand[1] ):
                                                lCommand[1] = string.replace( lCommand[1], "  ", " " )
                                        lVals = re.split( " ", lCommand[1] )
                                        self.lfTransformB = []
                                        for sVal in lVals:
                                                self.lfTransformB.append( float( sVal ) )
                                        #print(self.lfTransformB)
                                elif ( "basis function support" in lCommand[0].lower() ):
                                        lCommand[1] = lCommand[1].strip()
                                        if ( "global" in lCommand[1].lower() ):
                                                self.bGlobal = True
                                        elif ( "local" in lCommand[1].lower() ):
                                                self.bGlobal = False
                                        else:
                                                print(" ERROR: Basis Function Support should be either 'global' or 'local', instead of")
                                                print(lCommand[1])
                                                print("        at line {0:3d}".format( iCountLines ))
                                                return False
                                elif ( "points type" in lCommand[0].lower() ):
                                        self.sPointsType = lCommand[1].lower().strip()
                                        if ( not (self.sPointsType in lTasGridBuilderKnownPoints ) ):
                                                print(" ERROR: unrecognized points type")
                                                print(self.sPointsType)
                                                print("        at line {0:3d}".format( iCountLines ))
                                                print("")
                                                print("  acceptable points types are: ")
                                                for sS in lTasGridBuilderKnownPoints:
                                                        print("'{0:s}'".format( sS ) )
                                                print("")
                                                return False
                                elif ( "type parameter alpha" in lCommand[0].lower() ):
                                        self.fAlpha = float( lCommand[1] )
                                elif ( "type parameter beta" in lCommand[0].lower() ):
                                        self.fBeta = float( lCommand[1] )
                                elif ( "custom points file name" in lCommand[0].lower() ):
                                        self.sCustomFilename = lCommand[1]
                                elif ( "smolyak level" in lCommand[0].lower() ):
                                        self.iLevel = int( lCommand[1] )
                                        if ( self.iLevel < 1 ):
                                                print(" ERROR: Smolyak Level must be a positive integer")
                                                print(lCommand[1])
                                                print("        at line {0:3d}".format( iCountLines ))
                                                return False
                                elif ( "quadrature exactness" in lCommand[0].lower() ):
                                        self.iQExact = int( lCommand[1] )
                                        if ( self.iQExact < 0 ):
                                                print(" ERROR: Quadrature Exactness must be a non-negative integer")
                                                print(lCommand[1])
                                                print("        at line {0:3d}".format( iCountLines ))
                                                return False
                                elif ( "interpolation exactness" in lCommand[0].lower() ):
                                        self.iIExact = int( lCommand[1] )
                                        if ( self.iIExact < 0 ):
                                                print(" ERROR: Interpolation Exactness must be a non-negative integer")
                                                print(lCommand[1])
                                                print("        at line {0:3d}".format( iCountLines ))
                                                return False
                                elif ( "hyperbolic level" in lCommand[0].lower() ):
                                        self.iHyperbolic = int( lCommand[1] )
                                        if ( self.iHyperbolic < 0 ):
                                                print(" ERROR: Interpolation Exactness must be a non-negative integer")
                                                print(lCommand[1])
                                                print("        at line {0:3d}".format( iCountLines ))
                                                return False
                                elif ( "full tensor level" in lCommand[0].lower() ):
                                        lCommand[1] = string.replace( lCommand[1], "\t", " " )
                                        while ( "  " in lCommand[1] ):
                                                lCommand[1] = string.replace( lCommand[1], "  ", " " )
                                        lVals = re.split( " ", lCommand[1] )
                                        self.liTensor = []
                                        for sVal in lVals:
                                                self.liTensor.append( int( sVal ) )
                                        
                                        for iL in self.liTensor:
                                                if ( iL < 0 ):
                                                        print(" ERROR: tensor levels must be non-negative integers")
                                                        print(lCommand[1])
                                                        print("        at line {0:3d}".format( iCountLines ))
                                                        return False
                                elif ( "anisotropic weights" in lCommand[0].lower() ):
                                        lCommand[1] = string.replace( lCommand[1], "\t", " " )
                                        while ( "  " in lCommand[1] ):
                                                lCommand[1] = string.replace( lCommand[1], "  ", " " )
                                        lVals = re.split( " ", lCommand[1] )
                                        self.liAnisotropy = []
                                        for sVal in lVals:
                                                self.liAnisotropy.append( int( sVal ) )
                                        
                                        for iL in self.liAnisotropy:
                                                if ( iL < 1 ):
                                                        print(" ERROR: the anisotropic weights must be positive integers")
                                                        print(lCommand[1])
                                                        print("        at line {0:3d}".format( iCountLines ))
                                                        return False
                                elif ( "initial level" in lCommand[0].lower() ):
                                        self.iInitialLevel = int( lCommand[1] )
                                elif ( "basis functions type" in lCommand[0].lower() ):
                                        self.sFunctionType = lCommand[1].lower().strip()
                                        if ( not( self.sFunctionType in lTasGridBuilderKnownTypes ) ):
                                                print(" ERROR: unknown function type")
                                                print(lCommand[1])
                                                print("        at line {0:3d}".format( iCountLines ))
                                                print("")
                                                print("  acceptable basis types are: ")
                                                for sS in lTasGridBuilderKnownTypes:
                                                        print("'{0:s}'".format( sS ) )
                                                print("")
                                                return False
                                elif ( "basis order" in lCommand[0].lower() ):
                                        self.iBasisOrder = int( lCommand[1] )
                                        if ( self.iBasisOrder < -1 ):
                                                print(" ERROR: basis order should be non-negative or -1 to indicate automatic maximal order")
                                                print(lCommand[1])
                                                print("        at line {0:3d}".format( iCountLines ))
                                                return False
                                elif ( "enable refinement" in lCommand[0].lower() ):
                                        if ( "yes" in lCommand[1].lower() ):
                                                self.bRefinement = True
                                        elif ( "no" in lCommand[1].lower() ):
                                                self.bRefinement = False
                                        else:
                                                print(" ERROR: unrecognized 'enable refinement' option, say 'yes' or 'no'")
                                                print(lCommand[1])
                                                print("        at line {0:3d}".format( iCountLines ))
                                                return False
                                elif ( "maximum samples" in lCommand[0].lower() ):
                                        self.iMaxPoints = int( lCommand[1] )
                                        if ( self.iMaxPoints < 1 ):
                                                print(" ERROR: maximum number of points must be a positive integer")
                                                print(lCommand[1])
                                                print("        at line {0:3d}".format( iCountLines ))
                                                return False
                                elif ( "tolerance" in lCommand[0].lower() ):
                                        self.fTol = float( lCommand[1] )
                                        if ( self.fTol < 0.0 ):
                                                print(" ERROR: tolerance cannot be negative")
                                                print(lCommand[1])
                                                print("        at line {0:3d}".format( iCountLines ))
                                                return False
                                elif ( "refinement type" in lCommand[0].lower() ):
                                        self.sRefinementType = lCommand[1].lower().strip()
                                        if ( not( self.sRefinementType in lTasGridBuilderKnownRefinementTypes ) ):
                                                print(" ERROR: unrecognized 'Refinement Type'")
                                                print(lCommand[1])
                                                print("        at line {0:3d}".format( iCountLines ))
                                                print("")
                                                print("  acceptable refinement types are: ")
                                                for sS in lTasGridBuilderKnownRefinementTypes:
                                                        print("'{0:s}'".format( sS ) )
                                                print("")
                                                return False
                                else:
                                        print(" ERROR: unrecognized command")
                                        print(lCommand[0])
                                        print("        at line {0:3d}".format( iCountLines ))
                                        print("  Hint: commands are case insensitive, but should be spelled correctly with at least one white space between words")
                                        return False
                        iCountLines += 1
                
                bSane = self.checkSaneInput( bSilent )
                if ( bSane ):
                        self.preProcess()
                return bSane
                
        def checkSaneInput( self, bSilent ):
                bPass = True
                if ( not( 'sModelExecutableFile' in dir(self) ) ):
                        print(" ERROR: 'Model Executable File' not specified")
                        bPass = False
                if ( not( 'sModelInputFile' in dir(self) ) ):
                        print(" ERROR: 'Model Input File' not specified")
                        bPass = False
                if ( not( 'sModelOutputFile' in dir(self) ) ):
                        print(" ERROR: 'Model Output File' not specified")
                        bPass = False
                if ( not( 'sTasgridExec' in dir(self) ) ):
                        print(" ERROR: 'TasGrid Executable' not specified")
                        bPass = False
                if ( not( 'sProjectName' in dir(self) ) ):
                        print(" ERROR: 'Project Name' not specified")
                        bPass = False
                if ( not( 'iNumDimensions' in dir(self) ) ):
                        print(" ERROR: 'Number of Inputs' not specified")
                        bPass = False
                if ( not( 'iNumOutputs' in dir(self) ) ):
                        print(" ERROR: 'Number of Outputs' not specified")
                        bPass = False
                if ( ('lfTransformA' in dir(self) ) and ( 'iNumDimensions' in dir(self) ) ):
                        if ( (len( self.lfTransformA ) != 1) and (len( self.lfTransformA ) != self.iNumDimensions) ):
                                print(" ERROR: wrong number of entries in transform A, use either 1 entry or 1 entry per input")
                                bPass = False
                if ( ('lfTransformB' in dir(self) ) and ( 'iNumDimensions' in dir(self) ) ):
                        if ( (len( self.lfTransformB ) != 1) and (len( self.lfTransformB ) != self.iNumDimensions) ):
                                print(" ERROR: wrong number of entries in transform B, use either 1 entry or 1 entry per input")
                                bPass = False
                if ( ('lfTransformA' in dir(self) ) and ( not('lfTransformB' in dir(self) ) ) ):
                        print(" ERROR: transform A specified, but transform B is not specified, use both or neither")
                        bPass = False
                if ( ('lfTransformB' in dir(self) ) and ( not('lfTransformA' in dir(self) ) ) ):
                        print(" ERROR: transform B specified, but transform A is not specified, use both or neither")
                        bPass = False
                if ( not( 'bGlobal' in dir(self) ) ):
                        print(" ERROR: 'Basis Function Support' not specified, 'global' or 'local' must be specified")
                        bPass = False
                else: # check the local/global settings
                        if ( self.bGlobal ):
                                if ( not( 'sPointsType' in dir(self) ) ):
                                        print(" ERROR: 'Points Type' not specified")
                                        bPass = False
                                else:
                                        if ( (self.sPointsType in lTasGridBuilderRequireAlpha) and ( not( 'fAlpha' in dir(self) ) ) ):
                                                print(" ERROR: 'Type parameter Alpha' is not specified, but it is required by points type: '{0:s}'".format(self.sPointsType) )
                                                bPass = False
                                        elif ( ('fAlpha' in dir(self)) and (not bSilent) ):
                                                print(" WARNING: 'Type parameter Alpha' is specified, but it is ignored by points type '{0:s}'".format(self.sPointsType) )
                                        if ( (self.sPointsType in lTasGridBuilderRequireBeta) and ( not( 'fBeta' in dir(self) ) ) ):
                                                print(" ERROR: 'Type parameter Beta' is not specified, but it is required by points type '{0:s}'".format(self.sPointsType) )
                                                bPass = False
                                        elif ( ('fBeta' in dir(self)) and (not bSilent) ):
                                                print(" WARNING: 'Type parameter Beta' is specified, but it is ignored by points type '{0:s}'".format(self.sPointsType) )
                                        if ( ("custom-tabulated" in self.sPointsType) and ( not( 'sCustomFilename' in dir(self) ) ) ):
                                                print(" ERROR: 'Custom Points File Name' is not specified, but it is required by points type: '{0:s}'".format(self.sPointsType) )
                                                bPass = False
                                        elif ( ('sCustomFilename' in dir(self)) and (not bSilent) ):
                                                print(" WARNING: 'Custom Points File Name' is specified, but it is ignored by points type '{0:s}'".format(self.sPointsType) )
                                        if ( ("custom-tabulated" in self.sPointsType) and ( 'sCustomFilename' in dir(self) ) ):
                                                if ( not os.path.exists( self.sCustomFilename ) ):
                                                        print(" ERROR: 'Custom Points File Name' is set to '{0:s}', but this file does not exist".format( self.sCustomFilename ) )
                                                        bPass = False
                                                
                                iLCount = 0
                                if ( 'iLevel' in dir(self) ):
                                        iLCount += 1
                                if ( 'iQExact' in dir(self) ):
                                        iLCount += 1
                                if ( 'iIExact' in dir(self) ):
                                        iLCount += 1
                                if ( 'iHyperbolic' in dir(self) ):
                                        iLCount += 1
                                if ( 'liTensor' in dir(self) ):
                                        iLCount += 1
                                if ( iLCount == 0 ):
                                        print(" ERROR: must specify one of 'Smolyak Level', 'Quadrature Exactness', 'Interpolation Exactness', 'Hyperbolic Level', or 'Full Tensor Level'")
                                        bPass = False
                                if ( iLCount > 1 ):
                                        print(" ERROR: specified")
                                        if ( 'iLevel' in dir(self) ):
                                                print("Smolyak Level: {0:3d}".format( self.iLevel ) )
                                        if ( 'iQExact' in dir(self) ):
                                                print("Quadrature Exactness: {0:3d}".format( self.iQExact ) )
                                        if ( 'iIExact' in dir(self) ):
                                                print("Interpolation Exactness: {0:3d}".format( self.iIExact ) )
                                        if ( 'iHyperbolic' in dir(self) ):
                                                print("Hyperbolic Level: {0:3d}".format( self.iHyperbolic ) )
                                        if ( 'liTensor' in dir(self) ):
                                                sS = ""
                                                for iL in self.liTensor:
                                                        sS += " {0:2d}".format( iL )
                                                print("Full Tensor Level: {0:s}".format( sS ) )
                                        print(" however, only one can be used, hence only one should be specified")
                                        bPass = False
                                if ( ('liTensor' in dir(self)) and ('iNumDimensions' in dir(self) ) ):
                                        if ( (len( self.liTensor ) != 1) and ( len( self.liTensor ) != self.iNumDimensions ) ):
                                                print(" ERROR: wrong number of entries in 'Full Tensor Level', use either 1 entry or 1 entry per input")
                                                bPass = False
                                if ( ('liAnisotropy' in dir(self)) and ('iNumDimensions' in dir(self) ) ):
                                        if ( len(self.liAnisotropy) != self.iNumDimensions ):
                                                print(" ERROR: wrong number of 'Anisotropic Weights', requires 1 entry per input")
                                                bPass = False
                                if ( ('liTensor' in dir(self)) and ('liAnisotropy' in dir(self)) ):
                                        print(" ERROR: 'Full tensor Level' conflicts with 'Anisotropic Weights', cannot specify 'Anisotropic Weights' for full tensor grids")
                                        bPass = False
                        else: # local
                                if ( not( 'iInitialLevel' in dir(self) ) ):
                                        print(" ERROR: 'Initial Level' not specified")
                                        bPass = False
                                if ( not( 'sFunctionType' in dir(self) ) ):
                                        print(" ERROR: 'Basis Functions Type' not specified")
                                        bPass = False
                                if ( not( 'iBasisOrder' in dir(self) ) ):
                                        if ( not bSilent ):
                                                print(" WARNING: 'Basis Order' not specified, defaulting to linear")
                                else:
                                        if ( ( 'sFunctionType' in dir(self) ) and ( "wavelet" in self.sFunctionType ) ):
                                                if ( (self.iBasisOrder != 1) and (self.iBasisOrder != 3) ):
                                                        print(" ERROR: 'Basis Functions Type: wavelet' can only use 'Basis Order' of 1 or 3")
                                                        bPass = False
                # Refinement
                if ( not('bRefinement' in dir( self )) ):
                        if ( ('sPointsType' in dir(self)) and (self.sPointsType in lTasGridBuilderGlobalRefine) and (not bSilent ) ):
                                if ( not ('liTensor' in dir(self)) ):
                                        print(" WARNING: 'Enable Refinement' not specified, defaulting to 'no'")
                        self.bRefinement = False
                elif ( self.bRefinement and ('bGlobal' in dir(self)) ):
                        if ( self.bGlobal ):
                                if ( ('sPointsType' in dir(self)) and ( not( self.sPointsType in lTasGridBuilderGlobalRefine ) ) ):
                                        print("ERROR: 'Enable Refinement' is set to 'yes', however, no refinement is possible to points type: '{0:s}'".format(self.sPointsType) )
                                        bPass = False
                                if ( 'liTensor' in dir(self) ):
                                        print("ERROR: 'Enable Refinement' is set to 'yes', however, refinement is not possible for full tensor grids" )
                                        bPass = False
                                if ( ('sRefinementType' in dir(self)) and ( not bSilent) ):
                                        print(" WARNING: 'Refinement Type' is set, but global grids allow for only one type of refinement, hence ignoring 'Refinement Type'")
                        else:
                                if ( ( not('sRefinementType' in dir(self))) and ( not bSilent) ):
                                        print(" WARNING: 'Refinement Type' not specified, defaulting to 'FDS'")
                                        self.sRefinementType = "fds"
                        if ( not( 'fTol' in dir(self) ) ):
                                print(" ERROR: 'Enable Refinement' is set to 'yes', however, 'Tolerance' not specified")
                                bPass = False
                elif ( self.bRefinement ):
                        if ( not( 'fTol' in dir(self) ) ):
                                print(" ERROR: 'Enable Refinement' is set to 'yes', however, 'Tolerance' not specified")
                                bPass = False
                        #if ( not( 'sRefinementType' in dir(self) ) ):
                                #print(" WARNING: 'Enable Refinement' is set to 'yes', however, 'Refinement Type' not specified, defaulting to FDS")
                                #self.sRefinementType = "fds"
                                
                if ( not( 'iMaxPoints' in dir( self ) ) ):
                        self.iMaxPoints = -1
                        #if ( self.bRefinement

                return bPass
                
        def preProcess( self ):
                if ( ('lfTransformA' in dir(self) ) and (len( self.lfTransformA ) == 1) and ( self.iNumDimensions > 1 ) ):
                        fVal = self.lfTransformA[0]
                        for iI in range( self.iNumDimensions-1 ):
                                self.lfTransformA.append( fVal )
                if ( ('lfTransformB' in dir(self) ) and (len( self.lfTransformB ) == 1) and ( self.iNumDimensions > 1 ) ):
                        fVal = self.lfTransformB[0]
                        for iI in range( self.iNumDimensions-1 ):
                                self.lfTransformB.append( fVal )
                if ( ('liTensor' in dir(self) ) and (len( self.liTensor ) == 1) and ( self.iNumDimensions > 1 ) ):
                        fVal = self.liTensor[0]
                        for iI in range( self.iNumDimensions-1 ):
                                self.liTensor.append( fVal )
                if ( ('liAnisotropy' in dir(self) ) and (len( self.liAnisotropy ) == 1) and ( self.iNumDimensions > 1 ) ):
                        fVal = self.liAnisotropy[0]
                        for iI in range( self.iNumDimensions-1 ):
                                self.liAnisotropy.append( fVal )
                # set strings
                self.sDim = "{0:1d}".format( self.iNumDimensions )
                self.sOut = "{0:1d}".format( self.iNumOutputs )
                
                if ( self.bGlobal ):
                        self.sOnedim = self.sPointsType

                        if ( 'fAlpha' in dir( self ) ):
                                self.sAlpha = "{0:2.10e}".format( self.fAlpha )
                        if ( 'fBeta' in dir( self ) ):
                                self.sBeta = "{0:2.10e}".format( self.fBeta )

                        if ( 'iLevel' in dir( self ) ):
                                self.sType = 'level'
                                self.sDepth = "{0:1d}".format( self.iLevel )
                        elif ( 'iQExact' in dir( self ) ):
                                self.sType = 'qexact'
                                self.sDepth = "{0:1d}".format( self.iQExact )
                        elif ( 'iIExact' in dir( self ) ):
                                self.sType = 'iexact'
                                self.sDepth = "{0:1d}".format( self.iIExact )
                        elif ( 'iHyperbolic' in dir( self ) ):
                                self.sType = 'hyperbolic'
                                self.sDepth = "{0:1d}".format( self.iHyperbolic )
                        elif ( 'liTensor' in dir( self ) ):
                                self.sType = 'tensor'
                else:
                        if ( "wavelet" in self.sFunctionType ):
                                self.sOnedim = "local-wavelet"
                        elif ( "polynomial-zero-boundary" in self.sFunctionType ):
                                self.sOnedim = "local-polynomial-zero"
                        else:
                                self.sOnedim = "local-polynomial"
                                
                        self.sOrder = "{0:1d}".format( self.iBasisOrder )
                        
                if ( self.bRefinement ):
                        self.sTolerance = "{0:2.10e}".format( self.fTol )
                        
                self.sGridFile = self.sProjectName + "_MainGridFile.grid"
                self.sQFile = self.sProjectName + "_temp_QuadratureFile.matrix"
                self.sVFile = self.sProjectName + "_temp_ValuesFile.matrix"
                self.sWFile = self.sProjectName + "_temp_WeightsFile.matrix"
                self.sPFile = self.sProjectName + "_temp_PointsFile.matrix"
                self.sRFile = self.sProjectName + "_temp_ResultFile.matrix"
                self.sCFile = self.sProjectName + "_temp_CustomFile.matrix"
                
        def checkExistingProject( self ):
                bPassNames = True
                if ( os.path.exists( self.sGridFile ) ):
                        print(" ERROR: file '{0:s}' already exists, remove old files before starting a new project".format( self.sGridFile ) )
                        bPassNames = False
                if ( os.path.exists( self.sQFile ) ):
                        print(" ERROR: file '{0:s}' already exists, remove old files before starting a new project".format( self.sQFile ) )
                        bPassNames = False
                if ( os.path.exists( self.sVFile ) ):
                        print(" ERROR: file '{0:s}' already exists, remove old files before starting a new project".format( self.sVFile ) )
                        bPassNames = False
                if ( os.path.exists( self.sWFile ) ):
                        print(" ERROR: file '{0:s}' already exists, remove old files before starting a new project".format( self.sWFile ) )
                        bPassNames = False
                if ( os.path.exists( self.sPFile ) ):
                        print(" ERROR: file '{0:s}' already exists, remove old files before starting a new project".format( self.sPFile ) )
                        bPassNames = False
                if ( os.path.exists( self.sRFile ) ):
                        print(" ERROR: file '{0:s}' already exists, remove old files before starting a new project".format( self.sRFile ) )
                        bPassNames = False
                        
                return bPassNames
                        
        def deleteExistingProject( self ):
                if ( os.path.exists( self.sGridFile ) ):
                        os.system("rm -fr {0:s}".format( self.sGridFile ) )
                if ( os.path.exists( self.sQFile ) ):
                        os.system("rm -fr {0:s}".format( self.sQFile ) )
                if ( os.path.exists( self.sVFile ) ):
                        os.system("rm -fr {0:s}".format( self.sVFile ) )
                if ( os.path.exists( self.sWFile ) ):
                        os.system("rm -fr {0:s}".format( self.sWFile ) )
                if ( os.path.exists( self.sPFile ) ):
                        os.system("rm -fr {0:s}".format( self.sPFile ) )
                if ( os.path.exists( self.sRFile ) ):
                        os.system("rm -fr {0:s}".format( self.sRFile ) )
                
        def makeGrid( self ):
                sCommand  = self.sTasgridExec + " -makegrid"
                
                sCommand += " -gridfile " + self.sGridFile
                sCommand += " -outputfile " + self.sPFile
                
                sCommand += " -dimensions " + self.sDim
                sCommand += " -outputs " + self.sOut
                
                sCommand += " -onedim " + self.sOnedim
                
                if ( 'sDepth' in dir( self ) ):
                        sCommand += " -depth " + self.sDepth
                
                if ( self.bGlobal ):
                        sCommand += " -type " + self.sType
                        
                        if ( self.sPointsType in lTasGridBuilderRequireAlpha ):
                                sCommand += " -alpha " + self.sAlpha
                        if ( self.sPointsType in lTasGridBuilderRequireBeta ):
                                sCommand += " -beta " + self.sBeta
                                
                        # write anisotropy file
                        if ( 'liAnisotropy' in dir( self ) ):
                                lliMatrix = []
                                for iL in self.liAnisotropy:
                                        lL = [ iL ]
                                        lliMatrix.append( lL )
                                self.writeMatrix( lliMatrix, self.sWFile )
                                sCommand += " -anisotropyfile " + self.sWFile
                        
                        # tensor weights file
                        if ( "tensor" in self.sType ):
                                lliMatrix = []
                                for iL in self.liTensor:
                                        lL = [ iL ]
                                        lliMatrix.append( lL )
                                self.writeMatrix( lliMatrix, self.sWFile )
                                sCommand += " -anisotropyfile " + self.sWFile
                                
                        if ( "custom-rule" in self.sOnedim ):
                                sCommand += " -customrulefile " + self.sCustomFilename
                                
                        if ( "lfTransformA" in dir( self ) ):
                                lliMatrix = []
                                for iI in range( self.iNumDimensions ):
                                        lL = [ self.lfTransformA[iI], self.lfTransformB[iI] ]
                                        lliMatrix.append( lL )
                                self.writeMatrix( lliMatrix, self.sVFile )
                                sCommand += " -inputfile " + self.sVFile
                                
                else:
                        sCommand += " -order " + self.sOrder
                        
                
                if ( bDebugInfo ):
                        print(" Calling tasgrid with command")
                        print(sCommand)
                        
                os.system(sCommand)
                
                if ( self.bGlobal ):
                        if ( ("tensor" in self.sType) or ("liAnisotropy" in dir( self )) ):
                                os.system("rm -fr {0:s}".format( self.sWFile ) )
                                
                        if ( "lfTransformA" in dir( self ) ):
                                os.system("rm -fr {0:s}".format( self.sVFile ) )

        def callTheModel( self, bRecover, bSilent ):
                if ( os.path.exists( self.sPFile ) ):
                        llfPoints = self.readMatrix( self.sPFile )
                else: # points file is missing for some reason, get a new one from tasgrid
                        sCommand = self.sTasgridExec + " -getneededpoints"
                        sCommand += " -gridfile " + self.sGridFile
                        sCommand += " -outputfile " + self.sPFile
                        
                        os.system( sCommand )
                        
                        llfPoints = self.readMatrix( self.sPFile )
                
                llfVals = []
                if ( bRecover ):
                        if ( os.path.exists( self.sVFile ) ):
                                llfVals = self.readMatrix( self.sVFile )
                        elif ( not bSilent ):
                                print(" WARNING: cannot find a values file, maybe the crash happened at the first sample, restarting the last batch of samples from scratch")
                                
                while ( len( llfVals ) < len( llfPoints ) ):
                        iI = len( llfVals )
                        lfPoint = llfPoints[iI]
                        
                        self.writeInputFile( lfPoint, self.sModelInputFile )
                        
                        os.system( self.sModelExecutableFile )
                        
                        lfValue = self.readOutputFile( self.sModelOutputFile )
                        llfVals.append( lfValue )
                        
                        self.writeMatrix( llfVals, self.sVFile )
                        
                if ( len( llfVals ) > 0 ):
                        sCommand = self.sTasgridExec + " -loadvalues"
                        sCommand += " -gridfile " + self.sGridFile
                        sCommand += " -inputfile " + self.sVFile
                        
                        if ( bDebugInfo ):
                                print(" Calling tasgrid with command")
                                print(sCommand)
                                
                        os.system( sCommand )
                
                os.system("rm -fr {0:s}".format( self.sVFile ))
                os.system("rm -fr {0:s}".format( self.sPFile ))
                
                
        def setRefinement( self ):
                sCommand = self.sTasgridExec + " -refine"
                sCommand += " -gridfile " + self.sGridFile
                sCommand += " -outputfile " + self.sPFile
                sCommand += " -tolerance " + self.sTolerance
                
                if ( not self.bGlobal ):
                        sCommand += " -refinement " + self.sRefinementType
                        
                if ( bDebugInfo ):
                        print(" Calling tasgrid with command")
                        print(sCommand)
                        
                os.system( sCommand )
                
                if ( bDebugInfo ):
                        sCommand = self.sTasgridExec + " -s"
                        sCommand += " -gridfile " + self.sGridFile
                        os.system( sCommand )
                
                llfPoints = self.readMatrix( self.sPFile )
                
                if ( len( llfPoints ) > 0 ):
                        iNumCurrentPoints = self.getNumPoints()
                        if ( iNumCurrentPoints + len( llfPoints ) > self.iMaxPoints ):
                                print(" WARNING: tolerance has not been reached, however, the next refinement iteration will add more points than the specified maximum number of points")
                                print("          if this is unacceptable, update the maximum number of points and/or tolerance and restart with the --recover command")
                                
                                sCommand = self.sTasgridExec + " -cancelrefine "
                                sCommand += " -gridfile " + self.sGridFile
                                
                                if ( bDebugInfo ):
                                        print(" Calling tasgrid with command")
                                        print(sCommand)
                                
                                os.system( sCommand )
                                
                                os.system("rm -fr {0:s}".format( self.sPFile ) )
                                
                                return False
                        else:
                                return True
                else:
                        os.system("rm -fr {0:s}".format( self.sPFile ) )
                        return False
                
        def getNumPoints( self ):
                sCommand = self.sTasgridExec + " -summary"
                sCommand += " -gridfile " + self.sGridFile
                sCommand += " > " + self.sRFile
                
                if ( bDebugInfo ):
                        print(" Calling tasgrid with command")
                        print(sCommand)
                        
                os.system( sCommand )
                
                
                F = open( self.sRFile, 'r' )
                lsRLines = F.readlines()
                F.close()
                
                os.system("rm -fr {0:s}".format( self.sRFile ))
                
                for sLine in lsRLines:
                        if ( "number of points" in sLine ):
                                lsL = re.split( ":", sLine )
                                iNumPoints = int( lsL[1] )
                                return iNumPoints

def main( argv ):
        bCommandPass = True
        for iI in range( len( argv ) -1 ):
                sCom = argv[iI+1]
                if ( not ( sCom in lTasGridBuilderKnownCLICommands ) ):
                        print(" ERROR: unknown command  {0:s}".format( sCom ) )
                        bCommandPass = False
        
        if ( not bCommandPass ):
                return
        
        bSilent = False
        if ( ("-s" in argv) or ("--silent" in argv) ):
                bSilent = True
        
        bVerbose = False
        if ( ("-v" in argv) or ("--verbose" in argv) ):
                bVerbose = True
                
        bReset = False
        if ( "--reset" in argv ):
                bReset = True
                
        bRecover = False
        if ( "--recover" in argv ):
                bRecover = True
                
        if ( bRecover and bReset ):
                print(" ERROR: cannot simultaneously reset and recover")
                return
                
        
        print("------------------------------------------------------")
        TSW = TasGridWrapper()
        
        if ( not TSW.readDescription( bSilent ) ):
                print("\n ERROR: cannot initialize from ProblemDescription.in file")
                return
                
        if ( bReset ):
                TSW.deleteExistingProject()
        elif ( not bRecover ):
                if ( not( TSW.checkExistingProject() ) ):
                        print(" ERROR: files associated with this project already exist, use either the '--reset' or '--recover' commands")
                        return
                
        if ( bVerbose ):
                pass
        
        if ( not bRecover ):
                TSW.makeGrid()
        
        # set the initial level        
        TSW.callTheModel( bRecover, bSilent )
        
        # refinement
        bKeepRefining = True
        if ( TSW.bRefinement ):
                while ( bKeepRefining ):
                        bKeepRefining = TSW.setRefinement()
                        if ( bKeepRefining ):
                                TSW.callTheModel( False, bSilent )
                        #bKeepRefining = False
        
        
        #print("------------------------------------------------------")
        print("TasGrid file {0:s} generated.".format( TSW.sGridFile ) )
        print("------------------------------------------------------")
        
        if ( bDebugInfo ):
                sCommand = TSW.sTasgridExec + " -summary"
                sCommand += " -gridfile " + TSW.sGridFile
                os.system(sCommand)
                
main( sys.argv )
