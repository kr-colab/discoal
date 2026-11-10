#include "unity.h"
#include "discoal.h"
#include "discoalFunctions.h"
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <errno.h>
#include <time.h>

// Test fixtures
float *testData;
int testDataSize;

#ifndef TEST_RUNNER_MODE
void setUp(void) {
    trajectoryFd = -1;
    currentTrajectory = NULL;
    trajectoryFileSize = 0;
    testData = NULL;
    testDataSize = 0;
}

void tearDown(void) {
    if (trajectoryFd != -1) {
        if (currentTrajectory && currentTrajectory != MAP_FAILED) {
            munmap(currentTrajectory, trajectoryFileSize);
        }
        close(trajectoryFd);
        trajectoryFd = -1;
        currentTrajectory = NULL;
    }
    if (testData) {
        free(testData);
        testData = NULL;
    }
}
#endif

/* A test trajectory of numSteps floats in an anonymous file, as
 * proposeTrajectory() writes one: returns an open descriptor. */
int createTestTrajectoryFd(int numSteps) {
    FILE *fp = tmpfile();
    TEST_ASSERT_NOT_NULL(fp);
    free(testData);
    testData = (float *)malloc(numSteps * sizeof(float));
    TEST_ASSERT_NOT_NULL(testData);
    testDataSize = numSteps;
    for (int i = 0; i < numSteps; i++) {
        testData[i] = (float)i / numSteps;  // Values from 0.0 to ~1.0
    }
    TEST_ASSERT_EQUAL(numSteps, fwrite(testData, sizeof(float), numSteps, fp));
    TEST_ASSERT_EQUAL(0, fflush(fp));
    int fd = dup(fileno(fp));
    TEST_ASSERT_NOT_EQUAL(-1, fd);
    fclose(fp);
    return fd;
}

// Test ensureTrajectoryCapacity (deprecated function)
void test_ensureTrajectoryCapacity_within_limit(void) {
    // Should not crash for reasonable sizes
    ensureTrajectoryCapacity(1000);
    ensureTrajectoryCapacity(1000000);
    ensureTrajectoryCapacity(499999999);  // Just under limit
}

void test_ensureTrajectoryCapacity_exceeds_limit(void) {
    // This test would exit the program, so we can't test it directly
    // Just verify the limit is what we expect
    TEST_ASSERT_EQUAL(500000000, TRAJSTEPSTART);
}

// Test mmapAcceptedTrajectory
void test_mmapAcceptedTrajectory_basic(void) {
    int numSteps = 1000;
    mmapAcceptedTrajectory(createTestTrajectoryFd(numSteps), numSteps);
    
    // Verify globals are set correctly
    TEST_ASSERT_NOT_EQUAL(-1, trajectoryFd);
    TEST_ASSERT_NOT_NULL(currentTrajectory);
    TEST_ASSERT_NOT_EQUAL(MAP_FAILED, currentTrajectory);
    TEST_ASSERT_EQUAL(numSteps * sizeof(float), trajectoryFileSize);
    
    // Verify we can read the data
    for (int i = 0; i < numSteps; i++) {
        TEST_ASSERT_FLOAT_WITHIN(0.001, testData[i], currentTrajectory[i]);
    }
}

void test_mmapAcceptedTrajectory_replaces_previous(void) {
    // Create and map first trajectory
    int numSteps1 = 500;
    mmapAcceptedTrajectory(createTestTrajectoryFd(numSteps1), numSteps1);
    
    float *oldPtr = currentTrajectory;
    TEST_ASSERT_NOT_NULL(oldPtr);
    TEST_ASSERT_EQUAL(numSteps1 * sizeof(float), trajectoryFileSize);
    
    // Map a second trajectory with a different size
    int numSteps2 = 750;
    mmapAcceptedTrajectory(createTestTrajectoryFd(numSteps2), numSteps2);
    
    // Verify new mapping is active
    TEST_ASSERT_NOT_NULL(currentTrajectory);
    TEST_ASSERT_EQUAL(numSteps2 * sizeof(float), trajectoryFileSize);
    
    // If mmap returns same address (rare but possible), verify content is different
    if (oldPtr == currentTrajectory) {
        // Content should still be from new file
        TEST_ASSERT_FLOAT_WITHIN(0.001, (float)(numSteps2-1) / numSteps2, currentTrajectory[numSteps2-1]);
    }
}

void test_mmapAcceptedTrajectory_large_file(void) {
    // Test with a larger trajectory
    int numSteps = 100000;
    mmapAcceptedTrajectory(createTestTrajectoryFd(numSteps), numSteps);
    
    // Spot check some values
    TEST_ASSERT_FLOAT_WITHIN(0.001, 0.0f, currentTrajectory[0]);
    TEST_ASSERT_FLOAT_WITHIN(0.001, 0.5f, currentTrajectory[numSteps/2]);
    TEST_ASSERT_FLOAT_WITHIN(0.001, 0.99999f, currentTrajectory[numSteps-1]);
}

// Test mmap read-only protection
void test_mmapAcceptedTrajectory_read_only(void) {
    int numSteps = 100;
    mmapAcceptedTrajectory(createTestTrajectoryFd(numSteps), numSteps);
    
    // The mapping should be read-only (MAP_PRIVATE + PROT_READ)
    // We can read but not write
    float value = currentTrajectory[0];
    TEST_ASSERT_FLOAT_WITHIN(0.001, 0.0f, value);
    
    // Note: We can't easily test write protection without causing a segfault
}

/* The mapped trajectory's file has no name in the file system, so no
 * exit path can leave it behind. */
void test_mapped_trajectory_file_is_anonymous(void) {
    int numSteps = 500;
    mmapAcceptedTrajectory(createTestTrajectoryFd(numSteps), numSteps);
    struct stat st;
    TEST_ASSERT_EQUAL(0, fstat(trajectoryFd, &st));
    TEST_ASSERT_EQUAL(0, st.st_nlink);
    TEST_ASSERT_EQUAL(numSteps * sizeof(float), st.st_size);
}

/* A trajectory with no steps is not mapped (mmap cannot map zero bytes). */
void test_mmapAcceptedTrajectory_no_steps(void) {
    mmapAcceptedTrajectory(createTestTrajectoryFd(0), 0);
    TEST_ASSERT_NULL(currentTrajectory);
    TEST_ASSERT_NOT_EQUAL(-1, trajectoryFd);
}

#ifndef TEST_RUNNER_MODE
int main(void) {
    UNITY_BEGIN();
    
    RUN_TEST(test_ensureTrajectoryCapacity_within_limit);
    RUN_TEST(test_ensureTrajectoryCapacity_exceeds_limit);
    RUN_TEST(test_mmapAcceptedTrajectory_basic);
    RUN_TEST(test_mmapAcceptedTrajectory_replaces_previous);
    RUN_TEST(test_mmapAcceptedTrajectory_large_file);
    RUN_TEST(test_mmapAcceptedTrajectory_read_only);
    RUN_TEST(test_mapped_trajectory_file_is_anonymous);
    RUN_TEST(test_mmapAcceptedTrajectory_no_steps);
    
    return UNITY_END();
}
#endif