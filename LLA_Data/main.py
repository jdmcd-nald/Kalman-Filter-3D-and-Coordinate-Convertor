# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.
import re

def main():
    # Use a breakpoint in the code line below to debug your script.
    data = []
    with open("Data.kml","r") as file:
        for lines in file:
            current_line = lines.rstrip().lstrip()
            if current_line.find("latitude>") != -1:
                re.sub('<[^<]+>', "",current_line)
                print(re.sub('<[^<]+>', "",current_line))
                data.append(float(re.sub('<[^<]+>', "",current_line)))
            if current_line.find("longitude>") != -1:
                re.sub('<[^<]+>', "",current_line)
                print(re.sub('<[^<]+>', "",current_line))
                data.append(float(re.sub('<[^<]+>', "", current_line)))
            if current_line.find("altitude>") != -1:
                re.sub('<[^<]+>', "",current_line)
                print(re.sub('<[^<]+>', "",current_line))
                data.append(float(re.sub('<[^<]+>', "", current_line)))

    with open("LLA_Data.csv", "w") as file:
        time  = 10
        for x in range(int(len(data)/3)):
            file.write(str(time)+","+str(data[x*3+1])+","+str(data[x*3+0])+","+str(data[x*3+2])+"\n")
            time = time+10



# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
