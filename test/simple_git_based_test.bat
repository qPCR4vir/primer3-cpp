cd ../bin
echo off
echo Testing primer3_core.exe

echo
echo Test 1 : check '-output' and 'input_file'
primer3_core.exe  --default_version=1 --strict --output=../test/primer_check_output ../test/primer_check_input
echo Use git to see how the file primer_check_output changed
pause

echo .
echo Test 2 : check ' output' and ' input_file'
primer3_core.exe  --default_version=1 --st > ../test/primer_check_output < ../test/primer_check_input
echo Use git to see how the file primer_check_output changed
REM pause

echo .
echo Test 3 : check incorrect flag and '-error_file'
primer3_core.exe --flag --err=../test/cmd_test3_output < ../test/primer_check_input
echo Use git to see how the file cmd_test3_output changed
REM pause

echo .
echo Test 4 : check nonexistent input file and '2 error_file'
primer3_core.exe invalid_input 2> ../test/cmd_test4_output
echo Use git to see how the file cmd_test4_output changed
REM pause

echo .
echo Test 5 : check that io_version=3 fails
primer3_core.exe --io_version=3 2> ../test/cmd_test5_output
echo Use git to see how the file cmd_test4_output changed


echo .
echo Testing  ntdpal.exe

echo .
echo Test 1 :   # Test error handling on over-long input sequence:
ntdpal.exe ACGTGTTCGTCGTAAATAACATGCTATATT GACGTAGACAACCCTGTGTTTAGCCTGCGTTTTGTGCCATCCTAATGCTTTACTAGATCACTGAGCCACCTCCCAAGGACTACACCTAGCGGTATTTCGTACATTAACTAGGATCCTTTTCCACATGGACTACAATGTCTGCCGAGCATGCGATGGGGTACCGCGCCCGCGCACATACGCGCGCAGAGCTTTTGGAGGCATACCTACACCGGCGAGGGGCTGCGGTTTATTGACACTGAAACGGGATAACGAGTCGCTGAATTGAGCCAAAAATATGCAAGCGTCACAAATTGTGACAAAAATTTTAAAGGAAAAATTAGACCATTGATTCTGAAGTGGTGCGTATAGGACCAGTCGTGGCAATGAGACCGATTTGAGTAGCACTAGCTCAAACACTGTCTGGGTCGCCATCAAGGCCACAAGAACTTAAGCAGCCGTCACCCTATAGAAGGTTAAGCGACGGTTAGGGCTTCTGGCAACGAAAGTTGTCGGTTCGTCCTGTGCCAACGTGTGGCAAAGTCTACTATGATTCGATTGTTGACGTGTCGACAGGCTGTTTCGCTGGATACCCCACCTTGATAATTTTTCTCGTCGAACGCTAGCAGTTTTTTTTTCAACGGCCCGGAATCTGTAAGAGGCCGTTGCAGGAACGCGTGTGTATGTAAATGCCCACTACTTCTGTTATGTACCCAAATGGCGTGCGGCGTGGATGTATAGTGTCGACCCTCCATAATCGGGCGGACGGTCGTGGGGTATGTATGATCTTCGGCACTGATTCGCCTCGAGTCTATATGTTCTTAATCCAGACCTTCGGGGAAAGCCTACTTTCCATCCGTTGTCTAGCGTCATGCCAGTGACTACTGTTGTATTGTCTGGTTCCTAAGATAGCCATGGATTCCGGACATCGACGATGCACAAGAGCGTTAGCGCTGGTGTGCAACGCAACGTCGCGAAGGCTGGGTTACAGCGTGATCTCCTGGCTGCACCCAGATGCAGAGGGACATACCTACGATGAATAGGTGCGTCTGTTTATAAACGCCCAATCCTAGCAAAAATCACAACTAAGACAGTGTATGGAAGACCCACCAGTTGTGGGCGAATGGTCAGGTATACAAGATCGTGTCAAGACGGAACTTAAGCTTCTGTGCGCTCTCCATGCGAGCTGGTACGTCTGGACGGCGAGGTATGAGTGAATGACCATCCATGGCAACTTTCGTGTTCTACGACAGATACGAGCTCGACGGACGACCTGGTGACCAGTAGTATATGCGCGTCCGTCGGCCAGACTTTCCAAACGCCCTTTCAACGAGATACATGCGAACACGCTACAATTTCTCGTTCCGTCTAAAGTCGATACTCGCAAGCCCAGGCCCGTTACTACAACGCTGTTAATAGGATCAGAAGGGCCATAAGACTTTGGCAGCGGTAGCTAGGAAAGTGATGGTTGTGATGGCCCTAGTAAGGAGTCAGCCATCTACCCAACTATTTGAATGGGACCATAGCCAAGGGACCCAGCTGTTCCTTAGAAACCTGGTGACTCCCTTAGCCAATTGTGTAACTTCGTGCGTGCCAGTATTACACCTATAATCACAAGACCCCTTCAATACGAGTCCTGTGGCGTAGTGTTCCATCAAAACAATCAAGAACAGATTTCCGGTCCCCGTTGTGTTGGGATCTAGCGGACGTTGTCGGTAGATCAATAACGTAAATGCGAATCGAAGTTCTCTGGCCTAAAACAACTGCGCGCAGGGCCTCCGGTCATTGCATCTTTCTTGTCTCTCGTGAGGGCGTGATTCGTTTACCTGGAGCGAGCCGGGCACAAGAGCTATGGATTATTGGCTGGTGCAAAAACCATTCTAGCTACAATTATACTCGCGTGTCGACGATAAGAGTGAAATCACTGCGTAGGCAAACTGCCGGGTCACCAAGAGAGGCTGATACCGCGGTTCACCC  l  >../test/dpal_tmp 2>../test/dpal-err_tmp;
echo Use git to see how the file dpal.tmp changed. Beging with "Error: Sequence 2 longer than DPAL_MAX_ALIGN and alignment is requested\n" ?

echo .
echo Test 2 :   # 'Default implementations + alignment'
ntdpal.exe   > ../test/dpal_output < ../test/dpal_input 2>../test/dpal_output-err;
echo Use git to see how the file dpal_output changed.

echo .
echo Test 3 :   # 'Default implementations + NO alignment 1'
ntdpal.exe -s  > ../test/dpal_score_output < ../test/dpal_input 2>../test/dpal_score_output-err";
echo Use git to see how the file dpal_score_output changed.
pause

echo .
echo Test 4 :   # 'Default implementations + NO alignment 2'
ntdpal.exe -s  > ../test/dpal_long_score_output < ../test/dpal_long_input 2>../test/dpal_long_score_output-err";
echo Use git to see how the file dpal_long_score_output changed.


echo .
echo Test 5 :   # 'Force _dpal_generic'
ntdpal.exe -s  > ../test/dpal_score_output < ../test/dpal_input 2>&1";
echo Use git to see how the file dpal_score_output changed.

echo .
echo now exit...
pause


