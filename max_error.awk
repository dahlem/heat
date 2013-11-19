#!/bin/awk -f
# Print the maximum value in the forth column of a specified file.
# The fields in this file must be separated by commas.
BEGIN {
    FS = ",";
    MAX = 0;
}{
    if (NR > 1) {
	if ($4 > MAX) {
	    MAX = $4;
	}
    }
}
 END {
     print MAX
 }
