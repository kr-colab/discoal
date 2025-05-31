#ifndef UNITY_H
#define UNITY_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

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
            printf("FAIL: Expected non-NULL pointer at %s:%d\n", __FILE__, __LINE__); \
            exit(1); \
        } \
    } while(0)

#define TEST_ASSERT_EQUAL_DOUBLE(expected, actual) \
    do { \
        if ((expected) != (actual)) { \
            printf("FAIL: Expected %f but got %f at %s:%d\n", (double)(expected), (double)(actual), __FILE__, __LINE__); \
            exit(1); \
        } \
    } while(0)

// Test Runner
#define RUN_TEST(func) \
    do { \
        printf("Running " #func "...\n"); \
        setUp(); \
        func(); \
        tearDown(); \
        printf("PASS: " #func "\n"); \
    } while(0)

// Unity Test Runner Functions
#define UNITY_BEGIN() \
    printf("\n\nUnity Test Runner\n================\n\n")

#define UNITY_END() \
    printf("\n\nAll tests passed!\n"); \
    0

#ifdef __cplusplus
}
#endif

#endif /* UNITY_H */ 