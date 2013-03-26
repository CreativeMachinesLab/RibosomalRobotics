/* ---------------------------------------------------
   FILE:     constants.h
	AUTHOR:   Josh Bongard
	DATE:     March 21, 2001
	FUNCTION: This class contains all of the constants
			    for this project.
 -------------------------------------------------- */

#ifndef _CONSTANTS_H
#define _CONSTANTS_H

int			RANDOM_SEED			= 0;

int			TOTAL_OBJECTS			= 12;

int			AGE_GAP				= 1;

int                     AGE_CEILINGS[20]                = {
1*AGE_GAP,2*AGE_GAP,3*AGE_GAP,5*AGE_GAP,8*AGE_GAP,
13*AGE_GAP,21*AGE_GAP,34*AGE_GAP,55*AGE_GAP,89*AGE_GAP,
144*AGE_GAP,233*AGE_GAP,377*AGE_GAP,610*AGE_GAP,987*AGE_GAP,
1597*AGE_GAP,2584*AGE_GAP,4181*AGE_GAP,6765*AGE_GAP,10946*AGE_GAP};

int			NUM_SLOTS			= 2;

double			DEFAULT_MUTATION_RATE		= 0.01;
int			NUM_HIDDEN_NEURONS		= 0;

double			SUCCESSFUL_TARGET_OBJECT_PROXIMITY = 0.67;

double			CORE_DIST_MAX			= 0.4;
double			LEFT_SHOULDER_DIST_MAX		= 0.5;
double			RIGHT_SHOULDER_DIST_MAX		= LEFT_SHOULDER_DIST_MAX;

// ----------------------------------------------------------------
//                   General constants
// ----------------------------------------------------------------

int			EMPTY				= 0;
int			RUNNING				= 1;
int			FINISHED			= 2;

int			DATA_FILE_BUFFER		= 1;

char			TEMP_FILENAME_PREFIX_NEW[100]	= "tmp\Files";
char			TEMP_FILENAME[100] 		= "temp.dat";
char			TAUS_OUT_FILENAME[100]   	= "Taus.dat";
char			WEIGHTS_OUT_FILENAME[100]   	= "Weights.dat";
char			OMEGAS_OUT_FILENAME[100]   	= "Omegas.dat";
char			SENSORS_OUT_FILENAME[100]	= "SensorWeights.dat";

char			BODY_OUT_FILENAME[100] 		= "Body.dat";
char			SENSOR_IN_FILENAME[100] 	= "Sensors.dat";

int			EVALUATION_LENGTH		= 1000;

// ----------------------------------------------------------------
//                   Body constants
// ----------------------------------------------------------------

double			TARGET_RADIUS			= 0.35;

double			TARGET_LWH			= 1.0;

int			TARGET_DISTANCE			= 10;

double			SEGMENT_LENGTH			= TARGET_LWH*2.0;
double			SEGMENT_WIDTH			= TARGET_LWH;
double			SEGMENT_HEIGHT			= TARGET_LWH/4.0;
double			CORE_RADIUS			= SEGMENT_HEIGHT/4.0;
double			ARM_LENGTH			= SEGMENT_LENGTH-SEGMENT_WIDTH/2.0;

double			SPINE_VERTICAL_JOINT_RANGE	= 20.0;
double			SPINE_HORIZONTAL_JOINT_RANGE	= 50.0;

double			LEG_VERTICAL_JOINT_RANGE	= 10.0;
double			LEG_HORIZONTAL_JOINT_RANGE	= 50.0;

double			MOTOR_FORCE			= 20;

double			ARM_HEIGHT			= 1.0;

double			ARM_RADIUS			= 0.1;
double			UPPERARM_LENGTH			= 1.0;
double			UPPERARM_JOINT_RANGE		= 30;

double			SHOULDER_RADIUS			= 0.1;
double			SHOULDER_JOINT_RANGE		= 30;

double			LOWERARM_LENGTH			= 1.0;
double			LOWERARM_JOINT_RANGE		= 30;

double			WRIST_JOINT_RANGE		= 30;

double			HAND_RADIUS			= 0.25;

double			NUM_FINGERS			= 4.0;
double			INTER_FINGER_ANGLE		= 2*3.14159/NUM_FINGERS;
double			FINGER_RADIUS			= 0.075;
double			MIN_FINGER_RADIUS		= 0.2*FINGER_RADIUS;
double			MAX_FINGER_RADIUS		= 4.0*FINGER_RADIUS;

double			FINGER_LENGTH			= 0.3;
double			MIN_FINGER_LENGTH		= 0.01;
double			MAX_FINGER_LENGTH		= 0.6;

double			UPPER_FINGER_JOINT_RANGE	= 90;
double			FINGER_JOINT_RANGE		= 60;

/*
double			FINGER_MOTOR_FORCE		= 10;
double			ARM_MOTOR_FORCE			= 100;
double			MOTOR_SPEED			= 3;
*/
double			FINGER_MOTOR_FORCE		= 10;
double			ARM_MOTOR_FORCE			= 100;
double			MOTOR_SPEED			= 3;

// ----------------------------------------------------------------
//                   Brain constants
// ----------------------------------------------------------------

int			NUM_ORIENTATION_SENSORS		= 2;
int			NUM_RANGE_SENSORS		= 3;
int			NUM_TOUCH_SENSORS		= 4;
int			NUM_ANGLE_SENSORS		= 10;

int			NUM_SENSORS			= NUM_RANGE_SENSORS + 
							  NUM_TOUCH_SENSORS + 
							  NUM_ANGLE_SENSORS +
							  NUM_ORIENTATION_SENSORS;

int			NUM_NEURONS			= 10;

double			TAU_MIN				= 0.0001;
double			TAU_MAX				= 1.0;
double			WEIGHT_MIN			= -16;
double			WEIGHT_MAX			= 16;
double			OMEGA_MIN			= -4;
double			OMEGA_MAX			= 4;
double			SENSOR_MIN			= WEIGHT_MIN;
double			SENSOR_MAX			= WEIGHT_MAX;

#endif
