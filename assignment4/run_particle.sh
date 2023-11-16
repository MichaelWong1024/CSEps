#!/bin/bash
print_expected_result() {
  case $1 in
    0.000000)
      cat <<- 'EOF'
        particle id: (xpos, ypos) (xvel, yvel): last_update : (cwall, cpart)
        particle 0: (5.000000, 6.000000) (-0.090000, -0.380000): 0.000000 : (0, 0)
        particle 1: (3.000000, 2.000000) (0.150000, 0.100000): 0.000000 : (0, 0)
        particle 2: (6.000000, 1.000000) (-0.200000, 0.150000): 0.000000 : (0, 0)
        particle 3: (13.000000, 12.000000) (0.400000, -0.100000): 0.000000 : (0, 0)
        particle 4: (10.000000, 8.000000) (-0.150000, 0.150000): 0.000000 : (0, 0)
EOF
    ;;
    2.750000)
      cat <<- 'EOF'
        particle 0: (5.000000, 6.000000) (-0.090000, -0.380000): 0.000000 : (0, 0)
        particle 1: (3.000000, 2.000000) (0.150000, 0.100000): 0.000000 : (0, 0)
        particle 2: (6.000000, 1.000000) (-0.200000, 0.150000): 0.000000 : (0, 0)
        particle 3: (14.100000, 11.725000) (-0.400000, -0.100000): 2.750000 : (1, 0)
        particle 4: (10.000000, 8.000000) (-0.150000, 0.150000): 0.000000 : (0, 0)
EOF
    ;;
    3.404136)
      cat <<- 'EOF'
        particle 0: (5.000000, 6.000000) (-0.090000, -0.380000): 0.000000 : (0, 0)
        particle 1: (3.506120, 2.337414) (-0.200000, 0.150000): 3.374136 : (0, 1)
        particle 2: (5.325173, 1.506120) (0.150000, 0.100000): 3.374136 : (0, 1)
        particle 3: (14.100000, 11.725000) (-0.400000, -0.100000): 2.750000 : (1, 0)
        particle 4: (10.000000, 8.000000) (-0.150000, 0.150000): 0.000000 : (0, 0)        
EOF
    ;;
    5.202699)
      cat <<- 'EOF'
        particle 0: (4.537157, 4.045774) (-0.200000, 0.150000): 5.142699 : (0, 1)
        particle 1: (3.152408, 2.602698) (-0.090000, -0.380000): 5.142699 : (0, 2)
        particle 2: (5.325173, 1.506120) (0.150000, 0.100000): 3.374136 : (0, 1)
        particle 3: (14.100000, 11.725000) (-0.400000, -0.100000): 2.750000 : (1, 0)
        particle 4: (10.000000, 8.000000) (-0.150000, 0.150000): 0.000000 : (0, 0)        
EOF
    ;;
    9.400325)
      cat <<- 'EOF'
        particle 0: (4.537157, 4.045774) (-0.200000, 0.150000): 5.142699 : (0, 1)
        particle 1: (2.772821, 1.000000) (-0.090000, 0.380000): 9.360325 : (1, 2)
        particle 2: (5.325173, 1.506120) (0.150000, 0.100000): 3.374136 : (0, 1)
        particle 3: (14.100000, 11.725000) (-0.400000, -0.100000): 2.750000 : (1, 0)
        particle 4: (10.000000, 8.000000) (-0.150000, 0.150000): 0.000000 : (0, 0)        
EOF
    ;;
    13.307501)
      cat <<- 'EOF'
        particle 0: (4.537157, 4.045774) (-0.200000, 0.150000): 5.142699 : (0, 1)
        particle 1: (2.772821, 1.000000) (-0.090000, 0.380000): 9.360325 : (1, 2)
        particle 2: (5.325173, 1.506120) (0.150000, 0.100000): 3.374136 : (0, 1)
        particle 3: (9.889000, 10.672250) (-0.150000, 0.150000): 13.277501 : (1, 1)
        particle 4: (8.008375, 9.991625) (-0.400000, -0.100000): 13.277501 : (0, 1)        
EOF
    ;;
    16.702553)
      cat <<- 'EOF'
        particle 0: (2.231186, 5.775252) (-0.090000, 0.380000): 16.672553 : (0, 2)
        particle 1: (2.114721, 3.778646) (-0.200000, 0.150000): 16.672553 : (1, 3)
        particle 2: (5.325173, 1.506120) (0.150000, 0.100000): 3.374136 : (0, 1)
        particle 3: (9.889000, 10.672250) (-0.150000, 0.150000): 13.277501 : (1, 1)
        particle 4: (8.008375, 9.991625) (-0.400000, -0.100000): 13.277501 : (0, 1)        
EOF
    ;;
    20.000000)
      cat <<- 'EOF'
      particle id: (xpos, ypos) (xvel, yvel): last_update : (cwall, cpart)
      particle 0: (1.931716, 7.039682) (-0.090000, 0.380000): 20.000000 : (0, 2)
      particle 1: (1.449232, 4.277764) (-0.200000, 0.150000): 20.000000 : (1, 3)
      particle 2: (7.819052, 3.168707) (0.150000, 0.100000): 20.000000 : (0, 1)
      particle 3: (8.880625, 11.680625) (-0.150000, 0.150000): 20.000000 : (1, 1)
      particle 4: (5.319375, 9.319375) (-0.400000, -0.100000): 20.000000 : (0, 1)        

EOF
    ;;
    *)
      echo "No expected result defined for time $1."
    ;;
  esac
}

for time in 0.000000 2.750000 3.404136 5.202699 9.400325 13.307501 16.702553 20.000000
do
  echo "Running ./particle normal.txt for time ${time}"
  # Capture the actual output
  actual_output=$(./particle normal.txt "$time")
  echo "Actual Output:"
  echo "$actual_output"
  echo ""

  # Print the expected result
  echo "Expected Output:"
  print_expected_result "$time"
  echo ""
done