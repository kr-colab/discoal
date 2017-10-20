import sys, gzip
import numpy as np

def sortedFlankingPositionsByDistToTargSite(targetPos, flankingPositionsToExamine, desiredNumPositions, physLen):
    i=1
    sortedFlankingPositions = []

    while len(sortedFlankingPositions) < desiredNumPositions:
        lPos = targetPos-i
        rPos = targetPos+i
        if lPos >= 0 and lPos in flankingPositionsToExamine:
            sortedFlankingPositions.append(lPos)
        if rPos < physLen and rPos in flankingPositionsToExamine and len(sortedFlankingPositions) < desiredNumPositions:
            sortedFlankingPositions.append(rPos)
        i += 1

    return sortedFlankingPositions

def getNearestEmptyPositions(donorPos, snpCountAtPos, physLen):
    numColliders = snpCountAtPos[donorPos]-1

    freeSlots = {}
    for pos in snpCountAtPos:
        if snpCountAtPos[pos] == 0:
            freeSlots[pos] = 1
    assert len(freeSlots) >= numColliders

    return sortedFlankingPositionsByDistToTargSite(donorPos, freeSlots, numColliders, physLen)

def resolveCollision(donorPos, snpCountAtPos, physLen):
    for recipientPos in getNearestEmptyPositions(donorPos, snpCountAtPos, physLen):
        snpCountAtPos[recipientPos] += 1
        assert snpCountAtPos[recipientPos] == 1
        snpCountAtPos[donorPos] -= 1

def msPositionsToIntegerPositions(positions, physLen):
    assert physLen >= len(positions)

    snpCountAtPos = {}
    for i in range(physLen):
        snpCountAtPos[i] = 0
    for position in positions:
        intPos = int(physLen*position)
        if intPos == physLen:
            intPos = physLen-1
        snpCountAtPos[intPos] += 1

    collisions = {}
    for pos in snpCountAtPos:
        if snpCountAtPos[pos] > 1:
            collisions[pos] = 1

    midPos = physLen/2
    collisionPositions = []
    midHasCollision=0
    if midPos in collisions:
        collisionPositions.append(midPos)
        midHasCollision=1
    collisionPositions += sortedFlankingPositionsByDistToTargSite(midPos, collisions, len(collisions)-midHasCollision, physLen)
    for pos in collisionPositions:
        resolveCollision(pos, snpCountAtPos, physLen)

    assert max(snpCountAtPos.values()) == 1
    newPositions = [x for x in sorted(snpCountAtPos) if snpCountAtPos[x] > 0]
    assert newPositions[0] >= 0 and newPositions[-1] < physLen

    return newPositions

def main():
    usageStr="""usage:
python infiniteSitesToFinite.py msOutFileName physLen outFileName

This script takes simulation results in ms-style format (whose path is the first argument: msOutFileName) and converts the positions from continuous values into descrete chromosomal locations ranging from [0, physLen-1], where physLen is the second argument. The script then emits the positions of polymorphisms, with one line for each rep in the simulation file. Output is written to outFileName (third argument). If outFileName ends with .npz, a matrix of zeros and ones is saved in npz format, associated with the key 'm'. The length of each row in the matrix is equal to physLen, and positions that are polymorphic are denoted with 1, and those that are monomorphic are 0, with one row for each simulation replicate. If outFileName ends with .gz then for each simulation rep a tab-separated list of integer locations of polymorphisms is written on a separate line, and this file is compressed in gzip format. Otherwise, the output is the same as for .gz but is not compressed. The number of polymorphisms in each rep must be <= physLen, if you know what's good for you!
"""

    if len(sys.argv) != 4:
        sys.exit(usageStr)

    msOutFileName, physLen, outFileName = sys.argv[1:]
    physLen = int(physLen)
    if msOutFileName.endswith(".gz"):
        fopen = gzip.open
    else:
        fopen = open

    newPositionLists = []
    with fopen(msOutFileName) as msOutFile:
        for line in msOutFile:
            if line.startswith("positions:"):
                positions = [float(x) for x in line.lstrip("positions:").strip().split()]
                newPositions = msPositionsToIntegerPositions(positions, physLen)
                newPositionLists.append(newPositions)

    if outFileName.endswith(".npz"):
        m = np.zeros((len(newPositionLists), physLen),dtype="bool")
        for i in range(len(newPositionLists)):
            for j in newPositionLists[i]:
                m[i, j] = 1
        np.savez(outFileName, m=m)
    else:
        if outFileName.endswith(".gz"):
            fopen = gzip.open
        else:
            fopen = open
        with fopen(outFileName, "w") as outFile:
            for newPositions in newPositionLists:
                outFile.write("\t".join([str(x) for x in newPositions]) + "\n")

if __name__ == "__main__":
    main()
