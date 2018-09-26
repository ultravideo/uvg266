file1 = "G:/tools/Kvazaar-vvc/bin/x64-Release/out.log"
file2 = "G:/tools/VVCSoftware_VTM_git/bin/vs15/msvc-19.15/x86_64/release/all.vtmbmsstats"

file2_line = ""
file2_values = ["-1"]

def check_values(val1, val2):
    if check_int(val1[1]) == False:
        return True
    #print("Matching..")
    if int(val1[1]) == int(val2[2]) and val1[2].strip() == val2[4].strip():
        return True
    return False

def check_int(s):
    s = str(s).strip()
    #print("{0} {1} 2".format(s, s.isdigit()))

    return s.isdigit()

def RepresentsInt(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


with open(file1, "r") as fd1:
    with open(file2, "r") as fd2:
        for file1_line in fd1:
            file1_values = file1_line.strip().split()
            #print(file1_values[0])

            #print("{0} {1}".format(file1_values[0], check_int(file1_values[0])))
            if check_int(file1_values[0]) == False:
                continue
            #print("ok")
            while 1:
                if int(file1_values[0]) == int(file2_values[0]):
                    if check_values(file1_values, file2_values) != True:
                        print("Value mismatch at " + file1_values[0])
                        exit(1)
                    break
                if int(file1_values[0]) < int(file2_values[0]):
                    break
                file2_line = fd2.readline()
                file2_values = file2_line.split(" ")
                while check_int(file2_values[0]) == False:
                    file2_line = fd2.readline()
                    file2_values = file2_line.split(" ")


