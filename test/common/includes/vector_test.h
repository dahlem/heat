/* Copyright (C) 2008 Dominik Dahlem <Dominik.Dahlem@gmail.com>                */
/*                                                                             */
/* This file is free software; as a special exception the author gives         */
/* unlimited permission to copy and/or distribute it, with or without          */
/* modifications, as long as this notice is preserved.                         */
/*                                                                             */
/* This program is distributed in the hope that it will be useful, but         */
/* WITHOUT ANY WARRANTY, to the extent permitted by law; without even the      */
/* implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    */

#ifndef __VECTOR_TEST_H__
#define __VECTOR_TEST_H__


#include <CUnit/CUnit.h>


void registerVectorTests();



void testDotProduct();
void testDaxpy();
void testDnrm2();
void testAdd();
void testScale();


static CU_TestInfo test_vector[] = {
    { "testDotProduct", testDotProduct },
    { "testDaxpy", testDaxpy },
    { "testDnrm2", testDnrm2 },
    { "testAdd", testAdd },
    { "testScale", testScale },
    CU_TEST_INFO_NULL,
};

static CU_SuiteInfo vector_suites[] = {
    { "TestVector", NULL, NULL, test_vector },
    CU_SUITE_INFO_NULL,
};

#endif
