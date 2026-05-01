#include <math.h>
#include <stdlib.h>
#include "shapes.h"

double sizeAt(int popID, double t) {
    Shape *s = &popShape[popID];
    switch (s->type) {
        case SHAPE_CONSTANT:
            return s->anchor_value;
        default:
            return 0.0;  /* other shapes implemented in subsequent tasks */
    }
}

double migAt(int srcPopID, int dstPopID, double t) {
    (void)srcPopID; (void)dstPopID; (void)t;
    return 0.0;
}

double integratedHazardSize(int popID, double t0, double T, int k) {
    (void)popID; (void)t0; (void)T; (void)k;
    return 0.0;
}

double integratedHazardMig(int srcPopID, int dstPopID, double t0, double T, int k) {
    (void)srcPopID; (void)dstPopID; (void)t0; (void)T; (void)k;
    return 0.0;
}

double drawWaitingTimeSize(int popID, double t0, double xi, int k) {
    (void)popID; (void)t0; (void)xi; (void)k;
    return -1.0;
}

double drawWaitingTimeMig(int srcPopID, int dstPopID, double t0, double xi, int k) {
    (void)srcPopID; (void)dstPopID; (void)t0; (void)xi; (void)k;
    return -1.0;
}
