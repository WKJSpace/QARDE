#ifndef QARDE_TOOLS_HPP
#define QARDE_TOOLS_HPP

#include <ap_int.h>
#include <ap_fixed.h>

// Auto-bitwidth calculation for unsigned
// Base case: Number == 0 requires 1 bit (if needed to represent zero itself)
template <unsigned int Number, int Count = 0>
struct BITWIDTH_COUNT {
    static const int value = BITWIDTH_COUNT<(Number >> 1), Count + 1>::value;
};

template <int Count>
struct BITWIDTH_COUNT<0, Count> {
    static const int value = Count;
};

#define BIT_WIDTH(flags) (BITWIDTH_COUNT<(flags)>::value)
#define DATA_TYPE(flags) ap_uint<BITWIDTH_COUNT<(flags)>::value>

// Auto-bitwidth calculation for signed
template <int Number, int Count = 0>
struct SIGNED_BITWIDTH_COUNT {
    static const int value = SIGNED_BITWIDTH_COUNT<(Number >> 1), Count + 1>::value;
};

// Base case
template <int Count>
struct SIGNED_BITWIDTH_COUNT<0, Count> {
    static const int value = Count + 1; // Add 1 for the sign bit
};

// Macro for signed width
#define SIGNED_BIT_WIDTH(val) (SIGNED_BITWIDTH_COUNT<((val < 0) ? -(val) : val)>::value)
#define SIGNED_DATA_TYPE(val) ap_int<SIGNED_BITWIDTH_COUNT<((val < 0) ? -(val) : val)>::value>

// Bitonic sorting
template <typename T>
struct ValIdx {
    T val;
    int idx;
};

template<typename T>
struct TabCompEntry {
    T val;
    int x;
    int y;
};

// Bitonic sort for 2 elements
template <typename T>
void bitonic_sort2(ValIdx<T> a[2], ValIdx<T> b[2], bool sign) {
// Clang-format off
    #pragma HLS INLINE
    // Clang-format on
    if (sign) {
        if (a[0].val > a[1].val) {
            b[0] = a[1];
            b[1] = a[0];
        } else {
            b[0] = a[0];
            b[1] = a[1];
        }
    } else {
        if (a[0].val > a[1].val) {
            b[0] = a[0];
            b[1] = a[1];
        } else {
            b[0] = a[1];
            b[1] = a[0];
        }
    }
}

// Recursive bitonic sort template
template <typename T, int Number>
struct bitonic_sort_inst {
    static void sub_sort(ValIdx<T> a[Number], ValIdx<T> b[Number], bool sign) {
// Clang-format off
		#pragma HLS INLINE
		// Clang-format on
        ValIdx<T> temp1[Number / 2];
        ValIdx<T> temp2[Number / 2];
        ValIdx<T> temp3[Number / 2];
        ValIdx<T> temp4[Number / 2];
        ValIdx<T> temp5[Number];
        ValIdx<T> temp6[Number];
// Clang-format off
		#pragma HLS ARRAY_PARTITION variable = temp1 complete dim = 1
		#pragma HLS ARRAY_PARTITION variable = temp2 complete dim = 1
		#pragma HLS ARRAY_PARTITION variable = temp3 complete dim = 1
		#pragma HLS ARRAY_PARTITION variable = temp4 complete dim = 1
		#pragma HLS ARRAY_PARTITION variable = temp5 complete dim = 1
		#pragma HLS ARRAY_PARTITION variable = temp6 complete dim = 1
        // Clang-format on

        sort_loop0:
        for (int i = 0; i < Number / 2; i++) {
// Clang-format off
			#pragma HLS UNROLL
        	// Clang-format on
            temp1[i] = a[i];
            temp2[i] = a[i + Number / 2];
        }

        bitonic_sort_inst<T, Number / 2>::sub_sort(temp1, temp3, 1);
        bitonic_sort_inst<T, Number / 2>::sub_sort(temp2, temp4, 0);

    sort_loop1:
        for (int i = 0; i < Number / 2; i++) {
// Clang-format off
			#pragma HLS UNROLL
			// Clang-format on
            temp5[i] = temp3[i];
            temp5[i + Number / 2] = temp4[i];
        }

        bitonic_sort_inst<T, Number / 2>::sub_merge(temp5, temp6, sign);

    sort_loop2:
        for (int i = 0; i < Number; i++) {
// Clang-format off
			#pragma HLS UNROLL
			// Clang-format on
            b[i] = temp6[i];
        }
    }


    static void sub_merge(ValIdx<T> a[2 * Number], ValIdx<T> b[2 * Number], bool sign) {
// Clang-format off
        #pragma HLS INLINE
        // Clang-format on
    	int i;
    	ValIdx<T> temp00[Number][2];
    	ValIdx<T> temp01[Number][2];
    	ValIdx<T> temp1[Number];
    	ValIdx<T> temp2[Number];
    	ValIdx<T> temp3[Number];
    	ValIdx<T> temp4[Number];
//Clang-format off
		#pragma HLS ARRAY_PARTITION variable = temp00 complete
		#pragma HLS ARRAY_PARTITION variable = temp01 complete
		#pragma HLS ARRAY_PARTITION variable = temp1 complete dim = 1
		#pragma HLS ARRAY_PARTITION variable = temp2 complete dim = 1
		#pragma HLS ARRAY_PARTITION variable = temp3 complete dim = 1
		#pragma HLS ARRAY_PARTITION variable = temp4 complete dim = 1
    	// Clang-format on

		merge_loop0:
		for (i = 0; i < Number; i++) {
//Clang-format off
			#pragma HLS UNROLL
			// Clang-format on
			temp00[i][0] = a[i];
			temp00[i][1] = a[i + Number];
		}

		merge_loop1:
		for (i = 0; i < Number; i++) {
//Clang-format off
			#pragma HLS UNROLL
			// Clang-format on
			bitonic_sort2<T>(temp00[i], temp01[i], sign);
		}

		merge_loop2:
		for (i = 0; i < Number; i++) {
//Clang-format off
			#pragma HLS UNROLL
			// Clang-format on
			temp1[i] = temp01[i][0];
			temp2[i] = temp01[i][1];
		}

		bitonic_sort_inst<T, Number / 2>::sub_merge(temp1, temp3, sign);
		bitonic_sort_inst<T, Number / 2>::sub_merge(temp2, temp4, sign);

		merge_loop3:
		for (i = 0; i < Number; i++) {
//Clang-format off
			#pragma HLS UNROLL
			// Clang-format on
			b[i] = temp3[i];
			b[i + Number] = temp4[i];
		}
    }
};


// Specialization for base case
template <typename T>
struct bitonic_sort_inst<T, 2> {
    static void sub_sort(ValIdx<T> a[2], ValIdx<T> b[2], bool sign) {
// Clang-format off
        #pragma HLS INLINE
        // Clang-format on
        bitonic_sort2<T>(a, b, sign);
    }

    static void sub_merge(ValIdx<T> a[4], ValIdx<T> b[4], bool sign) {
// Clang-format off
        #pragma HLS INLINE
        // Clang-format on
        ValIdx<T> temp1[2], temp2[2], temp3[2], temp4[2];
        ValIdx<T> temp5[2], temp6[2], temp7[2], temp8[2];
// Clang-format off
        #pragma HLS ARRAY_PARTITION variable=temp1 complete
        #pragma HLS ARRAY_PARTITION variable=temp2 complete
        #pragma HLS ARRAY_PARTITION variable=temp3 complete
        #pragma HLS ARRAY_PARTITION variable=temp4 complete
        #pragma HLS ARRAY_PARTITION variable=temp5 complete
        #pragma HLS ARRAY_PARTITION variable=temp6 complete
        #pragma HLS ARRAY_PARTITION variable=temp7 complete
        #pragma HLS ARRAY_PARTITION variable=temp8 complete
        // Clang-format on

        temp1[0] = a[0], temp1[1] = a[2];
        temp2[0] = a[1], temp2[1] = a[3];

        bitonic_sort2<T>(temp1, temp3, sign);
        bitonic_sort2<T>(temp2, temp4, sign);

        temp5[0] = temp3[0], temp5[1] = temp4[0];
        temp6[0] = temp3[1], temp6[1] = temp4[1];

        bitonic_sort2<T>(temp5, temp7, sign);
        bitonic_sort2<T>(temp6, temp8, sign);

        b[0] = temp7[0], b[1] = temp7[1];
        b[2] = temp8[0], b[3] = temp8[1];
    }
};


// Top level function for finding the minimum and its index (True: Ascendent False: Descendant)
template <typename T, int Number>
void bitonicSort(T input[Number], ValIdx<T> sorted_validx[Number]) {
// Clang-format off
    #pragma HLS INLINE off
    // Clang-format on

    ValIdx<T> input_tmp[Number];
    ValIdx<T> output_tmp[Number];
// Clang-format off
    #pragma HLS ARRAY_PARTITION variable=input_tmp complete
    #pragma HLS ARRAY_PARTITION variable=output_tmp complete
    // Clang-format on

	for (int i = 0; i < Number; i++) {
// Clang-format off
		#pragma HLS UNROLL
		// Clang-format on
		input_tmp[i].val = input[i];
		input_tmp[i].idx = i;
	}

	bitonic_sort_inst<T, Number>::sub_sort(input_tmp, output_tmp, true);

	// Extract sorted values and indices
	for (int i = 0; i < Number; i++) {
// Clang-format off
		#pragma HLS UNROLL
		// Clang-format on
		sorted_validx[i].val = output_tmp[i].val;
		sorted_validx[i].idx = output_tmp[i].idx;
	}
}

// Find top N minimum values and their indices
template<typename T, int Number, int TopN>
void find_topN(const T input[Number], ValIdx<T> sorted_top[TopN], bool sign)
{
// Clang-format off
    #pragma HLS INLINE off
    // Clang-format on

    ValIdx<T> input_tmp[Number];
    ValIdx<T> output_tmp[Number];
// Clang-format off
    #pragma HLS ARRAY_PARTITION variable=input_tmp complete
    #pragma HLS ARRAY_PARTITION variable=output_tmp complete
    // Clang-format on

    // Initialize input_tmp with input values and indices
    for (int i = 0; i < Number; i++) {
// Clang-format off
        #pragma HLS UNROLL
        // Clang-format on
        input_tmp[i].val = input[i];
        input_tmp[i].idx = i;
    }

    bitonic_sort_inst<T, Number>::sub_sort(input_tmp, output_tmp, sign);

    // Extract top TopN values and their indices
    for (int i = 0; i < TopN; i++)
    {
// Clang-format off
        #pragma HLS UNROLL
        // Clang-format on
        sorted_top[i].val = output_tmp[i].val;
        sorted_top[i].idx = output_tmp[i].idx;
    }
}

template<typename T, int Number>
void insert_sort(ValIdx<T> updated_entry, int pos, ValIdx<T> sorted[Number])
{
// Clang-format off
    #pragma HLS INLINE off
	// CLang-format on

    ValIdx<T> tmp_entry = updated_entry;
    int i = pos;
    int i_final = i;

    // Phase 1: Original sequential logic for i_final
    int candidate_pos[Number - 1];
    bool candidate_valid[Number - 1];
// Clang-format off
    #pragma HLS ARRAY_PARTITION variable=candidate_pos complete
    #pragma HLS ARRAY_PARTITION variable=candidate_valid complete
    // Clang-format on

    // Phase 2: Evaluate all steps in parallel
    insert_sort_eval:
    for (int step = 1; step < Number; step++) {
// Clang-format off
        #pragma HLS UNROLL
    	// Clang-format on
        int left = i - step;
        int right = i + step;

        bool use_left  = (left >= 0) && (tmp_entry.val < sorted[left].val);
        bool use_right = (right < Number) && (tmp_entry.val > sorted[right].val);

        candidate_valid[step - 1] = use_left || use_right;
        candidate_pos[step - 1]   = use_left ? left : (use_right ? right : i);
    }

    // Phase 3: Sequentially select the last matched i_final
    i_final = i;  // default

    insert_sort_ifinal:
    for (int step = 1; step < Number; step++) {
// Clang-format off
        #pragma HLS UNROLL
    	// Clang-format on
        if (candidate_valid[step - 1]) {
            i_final = candidate_pos[step - 1];
        }
    }

    // Phase 4: In-place shift with a single loop
    ValIdx<T> original[Number];
// Clang-format off
	#pragma HLS ARRAY_PARTITION variable=original type=complete
    // Clang-format on

    // Copy original first to avoid RAW
    insert_sort_copy:
    for (int j = 0; j < Number; j++) {
// Clang-format off
        #pragma HLS UNROLL
    	// Clang-format on
        original[j] = sorted[j];
    }

    insert_sort_inplace:
    for (int j = 0; j < Number; j++) {
// Clang-format off
        #pragma HLS UNROLL
    	// Clang-format on
        if (i_final < i && j <= i && j > i_final) {
            sorted[j] = original[j - 1];  // right shift
        } else if (i_final > i && j >= i && j < i_final) {
            sorted[j] = original[j + 1];  // left shift
        } else if (j != i_final) {
            sorted[j] = original[j];     // retain
        }
    }

    // Phase 5: Insert
    sorted[i_final] = tmp_entry;
}

template<int Number>
int ap_countr_zero(ap_uint<Number> value) {
// Clang-format off
	#pragma HLS INLINE
	// Clang-format on
	ap_countr_zero:
    for (int i = 0; i < Number; i++) {
// Clang-format off
    	#pragma HLS UNROLL
    	// Clang-format on
        if (value[i]) return i;
    }
    return Number;  // No bits set
}

template <typename T, int Number>
int locateVal(const T arr[Number], T value) {
// Clang-format off
	#pragma HLS INLINE
	// Clang-format on
    int idx = -1;
    for (int i = 0; i < Number; i++) {
// Clang-format off
		#pragma HLS UNROLL
    	// Clang-format on
        if (arr[i] == value && idx == -1) {
            idx = i;
        }
    }
    return idx;
}

#endif // QARDE_TOOLS_HPP
