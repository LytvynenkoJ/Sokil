import numpy as np
import warnings
from constants import *

warnings.filterwarnings('ignore', category=Warning)

def hex_to_number(number):
	return np.uint64(int(number, base = 16))

def hex_array_to_number_array(number_array):
	return [hex_to_number(number) for number in number_array]

def number_to_hex(number):
	return hex(number)

def number_array_to_hex_array(number_array):
	return [number_to_hex(number) for number in number_array]

modulo = (1 << 64)

T = [hex_array_to_number_array(T_array) for T_array in T_RAW]

A_MUL = hex_array_to_number_array(A_MUL_RAW)
A_INV_MUL = hex_array_to_number_array(A_INV_MUL_RAW)

class modular_math:
	def get_byte_array(number, size):
		byte_array = np.empty(size)

		for i in range(size):
			byte_array[i] = (number >> np.uint64(64 * i)) & modulo

		return byte_array

	@staticmethod
	def inverse(number):
		return (number ^ np.uint64(modulo - 1)) % modulo
	
	@staticmethod
	def add(number_1, number_2):
		return (number_1 + number_2) % modulo

	@staticmethod
	def alpha_multiply(number):
		return (number << np.uint64(8)) ^ np.uint64(A_MUL[number >> np.uint64(56)])
		
	@staticmethod
	def alpha_inverse_multiply(number):
		return (number >> np.uint64(8)) ^ np.uint64(A_INV_MUL[number & np.uint64(255)])

	@staticmethod
	def non_linear_transform(number):
		result = np.uint64(0)

		for i in range(8):
			index = (number >> np.uint64(i * 8)) & np.uint64(255)
			
			result ^= np.uint64(T[i][index])
		
		return result

class Strumok:
	def _state_init_256(self, key, initial_vector):
		state = np.empty(16, dtype = np.uint64)

		state[15] = modular_math.inverse(key[0])
		state[14] = key[1]
		state[13] = modular_math.inverse(key[2])
		state[12] = key[3]
		state[11] = key[0]
		state[10] = modular_math.inverse(key[1])
		state[9] = key[2]
		state[8] = key[3]
		state[7] = modular_math.inverse(key[0])
		state[6] = modular_math.inverse(key[1])
		state[5] = key[2] ^ initial_vector[3]
		state[4] = key[3]
		state[3] = key[0] ^ initial_vector[2]
		state[2] = key[1] ^ initial_vector[1]
		state[1] = key[2]
		state[0] = key[3] ^ initial_vector[0]

		self.ctx["S"] = state

	def _state_init_512(self, key, initial_vector):
		state = np.empty(16, dtype = np.uint64)

		state[15] = key[0]
		state[14] = modular_math.inverse(key[1])
		state[13] = key[2]
		state[12] = key[3]
		state[11] = modular_math.inverse(key[7])
		state[10] = key[5]
		state[9] = modular_math.inverse(key[6])
		state[8] = key[4] ^ initial_vector[3]
		state[7] = modular_math.inverse(key[0])
		state[6] = key[1]
		state[5] = key[2] ^ initial_vector[2]
		state[4] = key[3]
		state[3] = key[4] ^ initial_vector[1]
		state[2] = key[5]
		state[1] = key[6]
		state[0] = key[7] ^ initial_vector[0]

		self.ctx["S"] = state

	def _next_init(self, i):
		self.ctx["S"][i] = \
			modular_math.alpha_multiply(self.ctx["S"][i]) ^ \
			modular_math.alpha_inverse_multiply(self.ctx["S"][(11 + i) % 16]) ^ \
			self.ctx["S"][(13 + i) % 16] ^ (
				(self.ctx["S"][(15 + i) % 16] + self.ctx["r"][0]) ^ self.ctx["r"][1]
		)

		new_r_0 = self.ctx["r"][1] + self.ctx["S"][(13 + i) % 16]
		self.ctx["r"][1] = modular_math.non_linear_transform(self.ctx["r"][0])
		self.ctx["r"][0] = new_r_0

	def _next(self, i):
		self.ctx["S"][i] = \
			modular_math.alpha_multiply(self.ctx["S"][i]) ^ \
			modular_math.alpha_inverse_multiply(self.ctx["S"][(11 + i) % 16]) ^ \
			self.ctx["S"][(13 + i) % 16]

		new_r_0 = self.ctx["r"][1] + self.ctx["S"][(13 + i) % 16]
		self.ctx["r"][1] = modular_math.non_linear_transform(self.ctx["r"][0])
		self.ctx["r"][0] = new_r_0

		return (self.ctx["S"][i] + self.ctx["r"][0]) ^ self.ctx["r"][1] ^ self.ctx["S"][(1 + i) % 16]

	def update_state(self):
		return np.array([self._next(i) for i in range(16)], dtype = np.uint64)

	def __init__(self, key, key_length, iv):
		self.ctx = {
			"key_size": key_length,
			"key": key,
			"iv": iv,
			"r": np.array([0, 0], dtype = np.uint64)
		}

		if (key_length == 256):
			self._state_init_256(key, iv)
		else:
			self._state_init_512(key, iv)

		for i in range(32):
			self._next_init(i % 16)
