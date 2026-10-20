/* (C) Copyright 2026- ECMWF. Licensed under the Apache Licence Version 2.0.
 * Streaming IEEE binary32 comparison: big-endian Fortran sequential spectra,
 * including continuation subrecords, or native little-endian diagnostic fields.
 */
#include <float.h>
#include <inttypes.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BLOCK (1024 * 1024)

static uint64_t values, different_words, different_values, failures, worst_index;
static double maximum, scaled_maximum, squared_sum;

static void fail(const char *message)
{
    fprintf(stderr, "%s\n", message);
    exit(1);
}

static uint32_t word(const unsigned char *bytes, int little)
{
    if (little)
        return (uint32_t)bytes[0] | (uint32_t)bytes[1] << 8 |
               (uint32_t)bytes[2] << 16 | (uint32_t)bytes[3] << 24;
    return (uint32_t)bytes[3] | (uint32_t)bytes[2] << 8 |
           (uint32_t)bytes[1] << 16 | (uint32_t)bytes[0] << 24;
}

static void compare(const unsigned char *left, const unsigned char *right, size_t bytes, int little)
{
    for (size_t index = 0; index < bytes; index += 4) {
        uint32_t a = word(left + index, little), b = word(right + index, little);
        if ((a & UINT32_C(0x7f800000)) == UINT32_C(0x7f800000) ||
            (b & UINT32_C(0x7f800000)) == UINT32_C(0x7f800000))
            fail("Nonfinite input: numerical acceptance cannot be established");
        if (a == b) {
            ++values;
            continue;
        }
        float ref, actual;
        memcpy(&ref, &a, sizeof(ref));
        memcpy(&actual, &b, sizeof(actual));
        double error = fabs((double)actual - ref);
        double scaled = error / fmax(1.0, fabs((double)ref));
        different_words += a != b;
        different_values += actual != ref;
        failures += error > 1.e-6 + 1.e-6 * fabs((double)ref);
        if (error > maximum) {
            maximum = error;
            worst_index = values;
        }
        scaled_maximum = fmax(scaled_maximum, scaled);
        squared_sum += error * error;
        ++values;
    }
}

static int marker(FILE *file, uint32_t *result)
{
    unsigned char bytes[4];
    size_t count = fread(bytes, 1, sizeof(bytes), file);
    if (!count && feof(file))
        return 0;
    if (count != sizeof(bytes))
        fail("Truncated Fortran record marker");
    *result = word(bytes, 0);
    return 1;
}

static uint64_t magnitude(uint32_t value)
{
    return value & UINT32_C(0x80000000) ? UINT64_C(0x100000000) - value : value;
}

int main(int argc, char **argv)
{
    int raw = argc == 4 && !strcmp(argv[1], "--raw-little");
    if (argc != 3 && !raw) {
        fprintf(stderr, "Usage: %s [--raw-little] REFERENCE ACTUAL\n", argv[0]);
        return 2;
    }
    if (sizeof(float) != 4 || FLT_RADIX != 2 || FLT_MANT_DIG != 24 ||
        FLT_MAX_EXP != 128 || FLT_MIN_EXP != -125)
        fail("IEEE binary32 host floats are required");
    FILE *left = fopen(argv[argc - 2], "rb"), *right = fopen(argv[argc - 1], "rb");
    if (!left || !right)
        fail("Cannot open both input files");
    unsigned char *a = malloc(BLOCK + 4), *b = malloc(BLOCK + 4);
    if (!a || !b)
        fail("Cannot allocate comparison buffers");
    uint64_t records = 0, subrecords = 0;
    size_t carry = 0;
    if (raw) {
        for (;;) {
            size_t na = fread(a + carry, 1, BLOCK, left);
            size_t nb = fread(b + carry, 1, BLOCK, right);
            if (ferror(left) || ferror(right))
                fail("Raw field read error");
            if (na != nb)
                fail("Raw field lengths differ");
            if (!na)
                break;
            size_t complete = (carry + na) / 4 * 4;
            compare(a, b, complete, 1);
            carry = carry + na - complete;
            memmove(a, a + complete, carry);
            memmove(b, b + complete, carry);
        }
    } else {
        for (;;) {
            uint32_t ma, mb;
            int ha = marker(left, &ma), hb = marker(right, &mb);
            if (ha != hb)
                fail("Fortran record counts differ");
            if (!ha)
                break;
            for (;;) {
                if (ma != mb || !magnitude(ma))
                    fail("Fortran record layouts differ or contain an empty record");
                uint64_t remaining = magnitude(ma), length = remaining;
                int more = !!(ma & UINT32_C(0x80000000));
                ++subrecords;
                while (remaining) {
                    size_t take = remaining > BLOCK ? BLOCK : (size_t)remaining;
                    if (fread(a + carry, 1, take, left) != take ||
                        fread(b + carry, 1, take, right) != take)
                        fail("Truncated Fortran record payload");
                    size_t complete = (carry + take) / 4 * 4;
                    compare(a, b, complete, 0);
                    carry = carry + take - complete;
                    memmove(a, a + complete, carry);
                    memmove(b, b + complete, carry);
                    remaining -= take;
                }
                uint32_t ta, tb;
                if (!marker(left, &ta) || !marker(right, &tb) || ta != tb || magnitude(ta) != length)
                    fail("Invalid Fortran trailing record marker");
                if (!more)
                    break;
                if (!marker(left, &ma) || !marker(right, &mb))
                    fail("Missing Fortran continuation subrecord");
            }
            if (carry)
                fail("Fortran logical record is not a sequence of binary32 values");
            ++records;
        }
    }
    if (carry || !values)
        fail("Input is empty or ends with a partial binary32 value");
    if (fclose(left) || fclose(right))
        fail("Cannot close input files");
    free(a);
    free(b);
    printf("{\"values\":%" PRIu64 ",\"bitwise_equal\":%s,\"differing_bit_patterns\":%" PRIu64
           ",\"differing_values\":%" PRIu64 ",\"tolerance_failures\":%" PRIu64
           ",\"max_absolute_error\":%.17g,\"max_scaled_error\":%.17g,\"rms_error\":%.17g"
           ",\"worst_flat_index\":%" PRIu64 ",\"records\":%" PRIu64 ",\"subrecords\":%" PRIu64 "}\n",
           values, different_words ? "false" : "true", different_words, different_values, failures,
           maximum, scaled_maximum, sqrt(squared_sum / values), worst_index, records, subrecords);
    return 0;
}
