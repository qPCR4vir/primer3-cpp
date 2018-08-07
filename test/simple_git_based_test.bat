cd ../bin
REM echo off
echo Testing primer3_core.exe

echo Test 1 : check '-output' and 'input_file'
primer3_core.exe  --default_version=1 --strict --output=../test/primer_check_output ../test/primer_check_input
echo Use git to see how the file primer_check_output changed
pause

echo Test 2 : check ' output' and ' input_file'
primer3_core.exe  --default_version=1 --st > ../test/primer_check_output < ../test/primer_check_input
echo Use git to see how the file primer_check_output changed
pause

echo Test 3 : check incorrect flag and '-error_file'
primer3_core.exe --flag --err=../test/cmd_test3_output < ../test/primer_check_input
echo Use git to see how the file cmd_test3_output changed
pause

echo Test 4 : check nonexistent input file and '2 error_file'
primer3_core.exe invalid_input 2> ../test/cmd_test4_output
echo Use git to see how the file cmd_test4_output changed
pause

echo Test 5 : check that io_version=3 fails
primer3_core.exe --io_version=3 2> ../test/cmd_test5_output
echo Use git to see how the file cmd_test4_output changed
pause
