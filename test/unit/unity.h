#ifndef UNITY_H
#define UNITY_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define UNITY_VERSION_MAJOR    2
#define UNITY_VERSION_MINOR    5
#define UNITY_VERSION_BUILD    1
#define UNITY_VERSION         ((UNITY_VERSION_MAJOR << 16) | (UNITY_VERSION_MINOR << 8) | UNITY_VERSION_BUILD)

#ifdef __cplusplus
extern "C"
{
#endif

// Test Suite Setup and Teardown
void setUp(void);
void tearDown(void);

// Test Assertions
#define TEST_ASSERT(condition) \
    do { \
        if (!(condition)) { \
            printf("FAIL: Assertion failed at %s:%d\n", __FILE__, __LINE__); \
            exit(1); \
        } \
    } while(0)

#define TEST_ASSERT_EQUAL(expected, actual) \
    do { \
        if ((expected) != (actual)) { \
            printf("FAIL: Expected %d but got %d at %s:%d\n", (int)(expected), (int)(actual), __FILE__, __LINE__); \
            exit(1); \
        } \
    } while(0)

#define TEST_ASSERT_FLOAT_EQUAL(expected, actual, delta) \
    do { \
        if (fabs((expected) - (actual)) > (delta)) { \
            printf("FAIL: Expected %f but got %f at %s:%d\n", (float)(expected), (float)(actual), __FILE__, __LINE__); \
            exit(1); \
        } \
    } while(0)

#define TEST_ASSERT_NULL(pointer) \
    do { \
        if ((pointer) != NULL) { \
            printf("FAIL: Expected NULL but got non-NULL at %s:%d\n", __FILE__, __LINE__); \
            exit(1); \
        } \
    } while(0)

#define TEST_ASSERT_NOT_NULL(pointer) \
    do { \
        if ((pointer) == NULL) { \
            printf("FAIL: Expected non-NULL but got NULL at %s:%d\n", __FILE__, __LINE__); \
            exit(1); \
        } \
    } while(0)

// Test Runner
#define RUN_TEST(func) \
    do { \
        printf("Running %s...\n", #func); \
        setUp(); \
        func(); \
        tearDown(); \
        printf("PASS: %s\n", #func); \
    } while(0)

#ifdef __cplusplus
}
#endif

#endif /* UNITY_H */ 