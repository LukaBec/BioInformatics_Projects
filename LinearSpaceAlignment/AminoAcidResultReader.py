def compare(v1,v2):
    start, end = 0, 0
    last = -1 #0 is match, 1 is mismatch/gap
    variant = ""
    for baseInd in range(len(v1)):
        if (v1[baseInd] == v2[baseInd]) and (last != 1):
            continue
        elif (v1[baseInd] != v2[baseInd]) and (last != 1):
            start = baseInd
            variant += v2[baseInd]
            last = 1
        elif (v1[baseInd] != v2[baseInd]) and (last == 1):
            variant += v2[baseInd]
        elif (v1[baseInd] == v2[baseInd]) and (last == 1):
            end = baseInd - 1
            last = 0
            print(variant, start, end)
            variant = ""


def read(file):
    with open(file, 'r') as f:
        all = [line.strip() for line in f.readlines()]
        return all[0], all[1]

def main():
    v1, v2 = "PLEASANTLY", "-MEA--N-LY" #read('Test.txt')
    print(v1)
    compare(v1,v2)

if __name__ == '__main__':
    main()
