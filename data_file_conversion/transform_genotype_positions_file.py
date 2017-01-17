#!/usr/bin/python3

class ConstructGenoTypeInformation ():

    # Init constructor.
    def __init__(self, snp_pos):

        # Set the location of the SNP positions file.
        self.snpPositions = snp_pos

    # Reads in the file, extracts information.
    def getGenoTypeInformation(self):

        # Opens the file.
        fh = open(self.snpPositions,"r")

        # Stores the content line by line in a variable.
        content = fh.readlines()

        # Close the file.
        fh.close()

        # A check to skip the initial header.
        skipFirstLine = 0

        # Handle for new file.
        newFile = self.constructGenoTypeInformationFile()

        # For every line do..
        for line in content:
            if skipFirstLine != 0:
                # Split the line on the tab ( note it's a csv file ).
                elements = line.split("\t")

                # Needed elements
                chr = elements[0]
                start = elements[2]
                end = int(elements[2]) +1
                id = elements[1]

                # Build new line.
                newLine = chr + "\t" + start + "\t" + str(end) + "\t" + id + "\n"

                self.writeToGenoTypeInformationFile(newFile, newLine)
            else:
                skipFirstLine += 1
                header = "chr\tstart\tend\tsnpId\n"
                self.writeToGenoTypeInformationFile(newFile, header)
        newFile.close()


    def constructGenoTypeInformationFile(self):
        return open("genoTypeInformation.tsv", "a")


    def writeToGenoTypeInformationFile(self, file, line):
        file.write(line)

def main():
    constructor = ConstructGenoTypeInformation(
        "/Users/molgenis/Dropbox/Erik Schutte Internship 2016/Data/eQTL-mapping-positions/CeD_43loci.txt")
    constructor.getGenoTypeInformation()
    constructor.constructGenoTypeInformationFile()
    
if "__main__" == __name__:
    main()
