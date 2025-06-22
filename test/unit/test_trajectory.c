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
char testTrajFilename[256];
float *testData;
int testDataSize;

#ifndef TEST_RUNNER_MODE
void setUp(void) {
    // Initialize globals
    trajectoryFd = -1;
    currentTrajectory = NULL;
    trajectoryFileSize = 0;
    testData = NULL;
    testDataSize = 0;
    
    // Create a unique test filename
    snprintf(testTrajFilename, sizeof(testTrajFilename), 
             "/tmp/test_traj_%d_%ld.tmp", getpid(), time(NULL));
}

void tearDown(void) {
    // Clean up any open trajectory
    if (trajectoryFd != -1) {
        if (currentTrajectory && currentTrajectory != MAP_FAILED) {
            munmap(currentTrajectory, trajectoryFileSize);
        }
        close(trajectoryFd);
        trajectoryFd = -1;
        currentTrajectory = NULL;
    }
    
    // Remove test file if it exists
    unlink(testTrajFilename);
    
    // Free test data
    if (testData) {
        free(testData);
        testData = NULL;
    }
}
#endif

// Helper function to create a test trajectory file
void createTestTrajectoryFile(const char *filename, int numSteps) {
    FILE *fp = fopen(filename, "wb");
    TEST_ASSERT_NOT_NULL(fp);
    
    testData = (float *)malloc(numSteps * sizeof(float));
    TEST_ASSERT_NOT_NULL(testData);
    testDataSize = numSteps;
    
    // Write some test data
    for (int i = 0; i < numSteps; i++) {
        testData[i] = (float)i / numSteps;  // Values from 0.0 to ~1.0
    }
    
    size_t written = fwrite(testData, sizeof(float), numSteps, fp);
    TEST_ASSERT_EQUAL(numSteps, written);
    
    fclose(fp);
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

// Test cleanupRejectedTrajectory
void test_cleanupRejectedTrajectory_existing_file(void) {
    // Create a test file
    FILE *fp = fopen(testTrajFilename, "w");
    TEST_ASSERT_NOT_NULL(fp);
    fprintf(fp, "test data\n");
    fclose(fp);
    
    // Verify file exists
    struct stat st;
    TEST_ASSERT_EQUAL(0, stat(testTrajFilename, &st));
    
    // Clean it up
    cleanupRejectedTrajectory(testTrajFilename);
    
    // Verify file is gone
    TEST_ASSERT_EQUAL(-1, stat(testTrajFilename, &st));
    TEST_ASSERT_EQUAL(ENOENT, errno);
}

void test_cleanupRejectedTrajectory_nonexistent_file(void) {
    // Should not crash when file doesn't exist
    cleanupRejectedTrajectory("/tmp/nonexistent_file_12345.tmp");
}

void test_cleanupRejectedTrajectory_null_filename(void) {
    // Should handle NULL gracefully
    cleanupRejectedTrajectory(NULL);
    cleanupRejectedTrajectory("");
}

// Test mmapAcceptedTrajectory
void test_mmapAcceptedTrajectory_basic(void) {
    int numSteps = 1000;
    createTestTrajectoryFile(testTrajFilename, numSteps);
    
    // Map the trajectory
    mmapAcceptedTrajectory(testTrajFilename, numSteps);
    
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
    createTestTrajectoryFile(testTrajFilename, numSteps1);
    mmapAcceptedTrajectory(testTrajFilename, numSteps1);
    
    float *oldPtr = currentTrajectory;
    TEST_ASSERT_NOT_NULL(oldPtr);
    TEST_ASSERT_EQUAL(numSteps1 * sizeof(float), trajectoryFileSize);
    
    // Create and map second trajectory with different size
    char secondFile[300];  // Larger buffer to avoid truncation warning
    snprintf(secondFile, sizeof(secondFile), "/tmp/test_traj2_%d_%ld.tmp", getpid(), time(NULL));
    int numSteps2 = 750;
    createTestTrajectoryFile(secondFile, numSteps2);
    mmapAcceptedTrajectory(secondFile, numSteps2);
    
    // Verify new mapping is active
    TEST_ASSERT_NOT_NULL(currentTrajectory);
    TEST_ASSERT_EQUAL(numSteps2 * sizeof(float), trajectoryFileSize);
    
    // If mmap returns same address (rare but possible), verify content is different
    if (oldPtr == currentTrajectory) {
        // Content should still be from new file
        TEST_ASSERT_FLOAT_WITHIN(0.001, (float)(numSteps2-1) / numSteps2, currentTrajectory[numSteps2-1]);
    }
    
    // Clean up second file
    unlink(secondFile);
}

void test_mmapAcceptedTrajectory_large_file(void) {
    // Test with a larger trajectory
    int numSteps = 100000;
    createTestTrajectoryFile(testTrajFilename, numSteps);
    
    mmapAcceptedTrajectory(testTrajFilename, numSteps);
    
    // Spot check some values
    TEST_ASSERT_FLOAT_WITHIN(0.001, 0.0f, currentTrajectory[0]);
    TEST_ASSERT_FLOAT_WITHIN(0.001, 0.5f, currentTrajectory[numSteps/2]);
    TEST_ASSERT_FLOAT_WITHIN(0.001, 0.99999f, currentTrajectory[numSteps-1]);
}

// Test mmap read-only protection
void test_mmapAcceptedTrajectory_read_only(void) {
    int numSteps = 100;
    createTestTrajectoryFile(testTrajFilename, numSteps);
    
    mmapAcceptedTrajectory(testTrajFilename, numSteps);
    
    // The mapping should be read-only (MAP_PRIVATE + PROT_READ)
    // We can read but not write
    float value = currentTrajectory[0];
    TEST_ASSERT_FLOAT_WITHIN(0.001, 0.0f, value);
    
    // Note: We can't easily test write protection without causing a segfault
}

// Test error handling (would require mocking or special setup)
void test_trajectory_filename_generation(void) {
    // Test that trajectory filenames are unique
    char name1[256], name2[256];
    
    snprintf(name1, sizeof(name1), "/tmp/discoal_traj_%d_%ld_%d.tmp", 
             getpid(), time(NULL), rand());
    
    usleep(1000);  // Small delay
    
    snprintf(name2, sizeof(name2), "/tmp/discoal_traj_%d_%ld_%d.tmp", 
             getpid(), time(NULL), rand());
    
    // Names should be different
    TEST_ASSERT_NOT_EQUAL(0, strcmp(name1, name2));
}

// Test cleanup behavior with multiple files
void test_multiple_trajectory_cleanup(void) {
    // Create multiple trajectory files
    char files[3][256];
    for (int i = 0; i < 3; i++) {
        // Ensure we don't overflow - limit base filename to leave room for suffix
        snprintf(files[i], sizeof(files[i]), "%.250s.%d", testTrajFilename, i);
        createTestTrajectoryFile(files[i], 100 + i * 50);
    }
    
    // Clean them up
    for (int i = 0; i < 3; i++) {
        cleanupRejectedTrajectory(files[i]);
    }
    
    // Verify all are gone
    struct stat st;
    for (int i = 0; i < 3; i++) {
        TEST_ASSERT_EQUAL(-1, stat(files[i], &st));
    }
}

// Test trajectory file persistence
void test_trajectory_file_persistence(void) {
    int numSteps = 500;
    createTestTrajectoryFile(testTrajFilename, numSteps);
    
    // Map the trajectory
    mmapAcceptedTrajectory(testTrajFilename, numSteps);
    
    // File should still exist on disk
    struct stat st;
    TEST_ASSERT_EQUAL(0, stat(testTrajFilename, &st));
    TEST_ASSERT_EQUAL(numSteps * sizeof(float), st.st_size);
    
    // Even after unmapping, file should persist (for accepted trajectories)
    if (currentTrajectory && currentTrajectory != MAP_FAILED) {
        munmap(currentTrajectory, trajectoryFileSize);
        currentTrajectory = NULL;
    }
    close(trajectoryFd);
    trajectoryFd = -1;
    
    TEST_ASSERT_EQUAL(0, stat(testTrajFilename, &st));
}

#ifndef TEST_RUNNER_MODE
int main(void) {
    UNITY_BEGIN();
    
    RUN_TEST(test_ensureTrajectoryCapacity_within_limit);
    RUN_TEST(test_ensureTrajectoryCapacity_exceeds_limit);
    RUN_TEST(test_cleanupRejectedTrajectory_existing_file);
    RUN_TEST(test_cleanupRejectedTrajectory_nonexistent_file);
    RUN_TEST(test_cleanupRejectedTrajectory_null_filename);
    RUN_TEST(test_mmapAcceptedTrajectory_basic);
    RUN_TEST(test_mmapAcceptedTrajectory_replaces_previous);
    RUN_TEST(test_mmapAcceptedTrajectory_large_file);
    RUN_TEST(test_mmapAcceptedTrajectory_read_only);
    RUN_TEST(test_trajectory_filename_generation);
    RUN_TEST(test_multiple_trajectory_cleanup);
    RUN_TEST(test_trajectory_file_persistence);
    
    return UNITY_END();
}
#endif