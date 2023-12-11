#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#define NORMAL 0
#define POS_INF 1
#define NEG_INF 2
#define NAN_VAL 3

double convert(char *bitstring, int totalBits, int mantissaBits);
void minmax(int totalBits, int mantissaBits);
void addHex(char* hex1, char* hex2, int totalBits, int mantissaBits);
int checkSpecialCase(char* binary, int totalBits, int mantissaBits);
void handleSpecialCases(int specialCase1, int specialCase2);
char* hexToBinary(char* hex, int totalBits);
double binaryToFloat(char* binary, int totalBits, int mantissaBits);
char* floatToBinary(double number, int totalBits, int mantissaBits);
char* binaryToHex(char* binary);

double convert(char *bitstring, int totalBits, int mantissaBits) {
    // Validate the length of the bitstring
    if ((int)strlen(bitstring) != totalBits) {
        printf("Error: Bitstring length does not match the specified total bits.\n");
        return 0.0;
    }

    // Extract the sign from the first bit of the bitstring
    int sign = (bitstring[0] == '0') ? 1 : -1;

    // Calculate the number of exponent bits and the bias
    int exponentBits = totalBits - mantissaBits - 1;
    int exponent = 0;
    int bias = pow(2, exponentBits - 1) - 1;

    // Extract and calculate the exponent value
    for (int i = 0; i < exponentBits; i++) {
        exponent = exponent * 2 + (bitstring[i + 1] - '0');
    }
    exponent -= bias;

    // Initialize mantissa value and calculate it by iterating over mantissa bits
    double mantissa = 1.0; // Start with 1.0 due to the implicit leading 1 in normalized numbers
    for (int i = 0; i < mantissaBits; i++) {
        mantissa += (bitstring[i + 1 + exponentBits] - '0') * pow(2, -1 - i);
    }

    // Check if all mantissa bits are zero (used for determining Infinity)
    bool isAllMantissaZero = true;
    for (int i = 1 + exponentBits; i < totalBits; i++) {
        if (bitstring[i] != '0') {
            isAllMantissaZero = false;
            break;
        }
    }
    
    // Determine if the number is Infinity or NaN
    bool isInfinity = exponent == (1 << exponentBits) - 1 && isAllMantissaZero;
    bool isNaN = exponent == (1 << exponentBits) - 1 && !isAllMantissaZero;
    if (isInfinity) {
        return sign * INFINITY; // Return Infinity with the correct sign
    }
    if (isNaN) {
        return NAN; // Return NaN
    }

    // Calculate and return the final floating-point value
    return sign * mantissa * pow(2, exponent);
}

void minmax(int totalBits, int mantissaBits) {
    int exponentBits = totalBits - mantissaBits - 1;
    char maxBitstring[totalBits + 1]; // +1 for null terminator
    char minBitstring[totalBits + 1];

    // Constructing bitstring for maximum value
    maxBitstring[0] = '0'; // Sign bit for positive number
    memset(maxBitstring + 1, '1', exponentBits - 1);
    maxBitstring[exponentBits] = '0'; // To avoid infinity
    memset(maxBitstring + exponentBits + 1, '1', mantissaBits);
    maxBitstring[totalBits] = '\0';

    // Constructing bitstring for minimum value
    minBitstring[0] = '0'; // Sign bit for positive number
    memset(minBitstring + 1, '0', exponentBits);
    memset(minBitstring + exponentBits + 1, '0', mantissaBits - 1);
    minBitstring[totalBits - 1] = '1'; // Smallest non-zero mantissa
    minBitstring[totalBits] = '\0';

    // Convert to base-10
    double minValue = convert(minBitstring, totalBits, mantissaBits);
    double maxValue = convert(maxBitstring, totalBits, mantissaBits);

    // Print the results
    printf("Min Positive: %.40f\n", minValue);
    printf("Max: %f\n", maxValue);
}

// void addHex(char* hex1, char* hex2, int totalBits, int mantissaBits) {
//     // Convert hexadecimal inputs to binary
//     char* binary1 = hexToBinary(hex1, totalBits);
//     char* binary2 = hexToBinary(hex2, totalBits);

//     // Check for special cases (+Inf, -Inf, NaN)
//     int specialCase1 = checkSpecialCase(binary1, totalBits, mantissaBits);
//     int specialCase2 = checkSpecialCase(binary2, totalBits, mantissaBits);

//     if (specialCase1 || specialCase2) {
//         handleSpecialCases(specialCase1, specialCase2);
//         free(binary1);
//         free(binary2);
//         return;
//     }

//     // Convert binary strings to floating-point numbers
//     double float1 = binaryToFloat(binary1, totalBits, mantissaBits);
//     double float2 = binaryToFloat(binary2, totalBits, mantissaBits);

//     // Perform addition
//     double sum = float1 + float2;

//     // Convert the result back to binary and hexadecimal
//     char* sumBinary = floatToBinary(sum, totalBits, mantissaBits);
//     char* sumHex = binaryToHex(sumBinary);

//     // Print results
//     printf("Binary 1: %s\n", binary1);
//     printf("Binary 2: %s\n", binary2);
//     printf("Value 1: %f\n", float1);
//     printf("Value 2: %f\n", float2);
//     printf("Sum: %f\n", sum);
//     printf("Binary Sum: %s\n", sumBinary);
//     printf("Hex Sum: %s\n", sumHex);

//     // Cleanup
//     free(binary1);
//     free(binary2);
//     free(sumBinary);
//     free(sumHex);
// }

void addHex(char* hex1, char* hex2, int totalBits, int mantissaBits) {
    // Convert hexadecimal inputs to binary
    char* binary1 = hexToBinary(hex1, totalBits);
    char* binary2 = hexToBinary(hex2, totalBits);

    // Check for special cases (like Infinity)
    int specialCase1 = checkSpecialCase(binary1, totalBits, mantissaBits);
    int specialCase2 = checkSpecialCase(binary2, totalBits, mantissaBits);

    // Print binary representations of inputs
    printf("Binary 1: %s\n", binary1);
    printf("Binary 2: %s\n", binary2);

    // Determine and print the values of each input
    double float1, float2;
    char* value1Str = (specialCase1 == POS_INF) ? "Inf" : ((specialCase1 == NEG_INF) ? "-Inf" : NULL);
    char* value2Str = (specialCase2 == POS_INF) ? "Inf" : ((specialCase2 == NEG_INF) ? "-Inf" : NULL);

    if (value1Str != NULL) {
        printf("Value 1: %s\n", value1Str);
        float1 = (specialCase1 == POS_INF) ? INFINITY : -INFINITY;
    } else {
        float1 = binaryToFloat(binary1, totalBits, mantissaBits);
        printf("Value 1: %f\n", float1);
    }

    if (value2Str != NULL) {
        printf("Value 2: %s\n", value2Str);
        float2 = (specialCase2 == POS_INF) ? INFINITY : -INFINITY;
    } else {
        float2 = binaryToFloat(binary2, totalBits, mantissaBits);
        printf("Value 2: %f\n", float2);
    }

    // Determine the sum and handle special case for Inf and -Inf
    double sum;
    char* sumBinary;
    char* sumHex;
    if ((specialCase1 == POS_INF && specialCase2 == NEG_INF) || (specialCase1 == NEG_INF && specialCase2 == POS_INF)) {
        sum = NAN; // Result is NaN
        sumBinary = floatToBinary(sum, totalBits, mantissaBits);
        sumHex = binaryToHex(sumBinary);
        printf("Sum: nan\n");
    } else {
        sum = float1 + float2;
        sumBinary = floatToBinary(sum, totalBits, mantissaBits);
        sumHex = binaryToHex(sumBinary);
        if (isinf(sum)) {
            printf("Sum: %s\n", (sum > 0) ? "Inf" : "-Inf");
        } else {
            printf("Sum: %f\n", sum);
        }
    }

    // Print binary and hex sums
    printf("Binary Sum: %s\n", sumBinary);
    printf("Hex Sum: %s\n", sumHex);

    // Cleanup
    free(binary1);
    free(binary2);
    free(sumBinary);
    free(sumHex);
}

int checkSpecialCase(char* binary, int totalBits, int mantissaBits) {
    int exponentBits = totalBits - mantissaBits - 1;
    char* exponent = (char*)malloc(exponentBits + 1);
    char* mantissa = (char*)malloc(mantissaBits + 1);

    if (!exponent || !mantissa) {
        // Handle memory allocation failure
        free(exponent);
        free(mantissa);
        return -1; // Indicates an error
    }

    strncpy(exponent, binary + 1, exponentBits); // Skip the sign bit
    exponent[exponentBits] = '\0';
    strncpy(mantissa, binary + 1 + exponentBits, mantissaBits);
    mantissa[mantissaBits] = '\0';

    // Check for all exponent bits set to 1 (indicative of special cases)
    int allOnes = 1;
    for (int i = 0; i < exponentBits; ++i) {
        if (exponent[i] != '1') {
            allOnes = 0;
            break;
        }
    }

    if (allOnes) {
        // Check mantissa to differentiate between Inf and NaN
        int allZeros = 1;
        for (int i = 0; i < mantissaBits; ++i) {
            if (mantissa[i] != '0') {
                allZeros = 0;
                break;
            }
        }

        free(exponent);
        free(mantissa);

        if (allZeros) {
            // +Inf or -Inf
            return binary[0] == '0' ? POS_INF : NEG_INF;
        } else {
            // NaN
            return NAN_VAL;
        }
    }

    free(exponent);
    free(mantissa);
    return NORMAL; // Not a special case
}

void handleSpecialCases(int specialCase1, int specialCase2) {
    if (specialCase1 == NAN_VAL || specialCase2 == NAN_VAL) {
        printf("Result: nan\n");
    } else if (specialCase1 == POS_INF || specialCase2 == POS_INF) {
        printf("Result: +Infinity\n");
    } else if (specialCase1 == NEG_INF || specialCase2 == NEG_INF) {
        printf("Result: -Infinity\n");
    } else {
        printf("Result: Undefined behavior\n");
    }
}

char* hexToBinary(char* hex, int totalBits) {
    // Allocate memory for the binary string
    char* binary = (char*)malloc(totalBits + 1);
    if (!binary) {
        printf("Memory allocation failed.\n");
        return NULL;
    }

    binary[0] = '\0'; // Start with an empty string

    // Iterate over each hex character and append its binary representation
    for (int i = 0; hex[i] != '\0'; i++) {
        switch (hex[i]) {
            case '0': strcat(binary, "0000"); break;
            case '1': strcat(binary, "0001"); break;
            case '2': strcat(binary, "0010"); break;
            case '3': strcat(binary, "0011"); break;
            case '4': strcat(binary, "0100"); break;
            case '5': strcat(binary, "0101"); break;
            case '6': strcat(binary, "0110"); break;
            case '7': strcat(binary, "0111"); break;
            case '8': strcat(binary, "1000"); break;
            case '9': strcat(binary, "1001"); break;
            case 'A': case 'a': strcat(binary, "1010"); break;
            case 'B': case 'b': strcat(binary, "1011"); break;
            case 'C': case 'c': strcat(binary, "1100"); break;
            case 'D': case 'd': strcat(binary, "1101"); break;
            case 'E': case 'e': strcat(binary, "1110"); break;
            case 'F': case 'f': strcat(binary, "1111"); break;
            default: 
                printf("Invalid hexadecimal character.\n");
                free(binary);
                return NULL;
        }
    }

    int currentLength = strlen(binary);
    if (currentLength > totalBits) {
        printf("Error: Binary representation exceeds total bits.\n");
        free(binary);
        return NULL;
    }

    // Sign extension
    char signBit = binary[0]; // Determine the sign bit from the first bit of the binary string
    while (currentLength < totalBits) {
        // Prepend the sign bit
        memmove(binary + 1, binary, currentLength + 1);
        binary[0] = signBit;
        currentLength++;
    }

    binary[totalBits] = '\0'; // Null-terminate the string
    return binary;
}

double binaryToFloat(char* binary, int totalBits, int mantissaBits) {
    if (!binary) return NAN; // Check for null pointer

    // Extract the sign from the first bit of the binary string
    int sign = (binary[0] == '0') ? 1 : -1;

    // Calculate the number of bits for the exponent
    int exponentBits = totalBits - mantissaBits - 1;
    int exponent = 0;
    int bias = pow(2, exponentBits - 1) - 1;

    // Extract and calculate the exponent value
    for (int i = 0; i < exponentBits; i++) {
        exponent = exponent * 2 + (binary[i + 1] - '0');
    }
    exponent -= bias;

    // Initialize mantissa value and calculate it by iterating over mantissa bits
    double mantissa = 1.0; // Start with 1.0 due to the implicit leading 1 in normalized numbers
    for (int i = 0; i < mantissaBits; i++) {
        mantissa += (binary[i + 1 + exponentBits] - '0') * pow(2, -1 - i);
    }

    // Check if all mantissa bits are zero (used for determining Infinity)
    bool isAllMantissaZero = true;
    for (int i = 1 + exponentBits; i < totalBits; i++) {
        if (binary[i] != '0') {
            isAllMantissaZero = false;
            break;
        }
    }
    
    // Determine if the number is Infinity or NaN
    bool isInfinity = exponent == (1 << exponentBits) - 1 && isAllMantissaZero;
    bool isNaN = exponent == (1 << exponentBits) - 1 && !isAllMantissaZero;
    if (isInfinity) {
        return sign * INFINITY; // Return Infinity with the correct sign
    }
    if (isNaN) {
        return NAN; // Return NaN
    }

    // Calculate and return the final floating-point value
    return sign * mantissa * pow(2, exponent);
}

char* floatToBinary(double number, int totalBits, int mantissaBits) {
    int exponentBits = totalBits - mantissaBits - 1;
    int bias = pow(2, exponentBits - 1) - 1;
    char* binary = (char*)malloc(totalBits + 1);

    if (!binary) {
        printf("Memory allocation failed.\n");
        return NULL;
    }

    // Handling sign bit
    binary[0] = (number < 0) ? '1' : '0';
    number = fabs(number);

    // Handling special cases: zero, infinity, and NaN
    if (number == 0) {
        // Set all bits to 0 for zero
        for (int i = 1; i < totalBits; i++) {
            binary[i] = '0';
        }
    } else if (isinf(number)) {
        // Set exponent bits to 1 and mantissa to 0 for infinity
        for (int i = 1; i <= exponentBits; i++) {
            binary[i] = '1';
        }
        for (int i = 1 + exponentBits; i < totalBits; i++) {
            binary[i] = '0';
        }
    } else if (isnan(number)) {
        // Set exponent bits to 1 and at least one mantissa bit to 1 for NaN
        for (int i = 1; i <= exponentBits; i++) {
            binary[i] = '1';
        }
        binary[1 + exponentBits] = '1';
        for (int i = 2 + exponentBits; i < totalBits; i++) {
            binary[i] = '0';
        }
    } else {
        // Normal case: Normalize the number and calculate exponent and mantissa

        // Normalize the number to get the exponent
        int exponent = 0;
        while (number >= 2.0) {
            number /= 2.0;
            exponent++;
        }
        while (number < 1.0) {
            number *= 2.0;
            exponent--;
        }

        // Adjust exponent with bias
        exponent += bias;

        // Convert exponent to binary
        for (int i = exponentBits; i > 0; i--) {
            binary[i] = (exponent % 2) ? '1' : '0';
            exponent /= 2;
        }

        // Convert mantissa to binary
        number -= 1.0; // Remove leading 1
        for (int i = 1 + exponentBits; i < totalBits; i++) {
            number *= 2.0;
            if (number >= 1.0) {
                binary[i] = '1';
                number -= 1.0;
            } else {
                binary[i] = '0';
            }
        }
    }

    binary[totalBits] = '\0'; // Null-terminate the string
    return binary;
}

char* binaryToHex(char* binary) {
    // Calculate the length of the binary string and the corresponding hex string length
    int binaryLength = strlen(binary);
    int hexLength = (binaryLength + 3) / 4; // Each hex digit represents 4 binary digits

    // Allocate memory for the hex string (+1 for the null-terminator)
    char* hex = (char*)malloc(hexLength + 1);
    if (!hex) {
        printf("Memory allocation failed.\n");
        return NULL;
    }

    hex[hexLength] = '\0'; // Null-terminate the string

    // Iterate over each group of 4 binary digits to convert to a hex digit
    for (int i = 0; i < hexLength; i++) {
        int value = 0; // Variable to store the decimal value of each 4-bit binary group
        int pow = 1; // Power of 2 (1, 2, 4, 8)
        for (int j = 0; j < 4; j++) {
            // Calculate the index of the binary digit
            int idx = binaryLength - (i * 4 + j) - 1;

            // Convert and add the binary digit value if the index is within the binary string
            if (idx >= 0) {
                value += (binary[idx] - '0') * pow;
            }
            pow *= 2; // Increment the power of 2 for the next bit
        }

        // Convert the decimal value to a hexadecimal digit and store in the hex string
        if (value < 10) {
            hex[hexLength - 1 - i] = '0' + value; // For 0-9
        } else {
            hex[hexLength - 1 - i] = 'A' + (value - 10); // For A-F
        }
    }

    return hex; // Return the hexadecimal string
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <totalBits> <mantissaBits> <operation> [additional arguments]\n", argv[0]);
        return 1;
    }

    int totalBits = atoi(argv[1]);
    int mantissaBits = atoi(argv[2]);

    if (totalBits < 8 || totalBits > 64 || mantissaBits < 1 || mantissaBits >= totalBits) {
        fprintf(stderr, "Invalid number of bits specified.\n");
        return 1;
    }

    if (strcmp(argv[3], "convert") == 0) {
        if (argc != 5) {
            fprintf(stderr, "Usage: %s <totalBits> <mantissaBits> convert <bitstring>\n", argv[0]);
            return 1;
        }
        double result = convert(argv[4], totalBits, mantissaBits);
        printf("%.10f\n", result);
    }
    else if (strcmp(argv[3], "minmax") == 0) {
        if (argc != 4) {
            fprintf(stderr, "Usage: %s <totalBits> <mantissaBits> minmax\n", argv[0]);
            return 1;
        }
        minmax(totalBits, mantissaBits);
    }
    else if (strcmp(argv[3], "addhex") == 0) {
        if (argc != 6) {
            fprintf(stderr, "Usage: %s <totalBits> <mantissaBits> addhex <hex1> <hex2>\n", argv[0]);
            return 1;
        }
        addHex(argv[4], argv[5], totalBits, mantissaBits);
    }
    else {
        fprintf(stderr, "Invalid operation specified.\n");
        return 1;
    }
}
