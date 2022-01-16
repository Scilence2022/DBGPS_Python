from DNAdroplet import DNADroplet
from fountain import Fountain
from utils import xor, randChunkNums
import random


class DNAFountain(Fountain):
	def __init__(self, data=b'', chunk_size=35, index=1, seed=1, droplet_head_length=16, droplet_tail_length=16, droplet_crc_length=8):
		self.droplet_head_length = droplet_head_length
		self.droplet_tail_length = droplet_tail_length
		self.droplet_crc_length = droplet_crc_length
		super(DNAFountain, self).__init__(data, chunk_size, index, seed)


	def DNAdroplet(self):
		self.updateIndex()
		aDroplet = DNADroplet(b'', self.index, self.num_of_chunks, self.droplet_head_length, self.droplet_tail_length, self.seed, self.droplet_crc_length)
		aDroplet.data_len = self.chunk_size * 4

		#aDroplet.degree = self.degrees[self.index % len(self.degrees)]
		aDroplet.degree = self.degrees[self.index -1]

		chunk_nums = aDroplet.get_chunk_nums()
		data = None
		for num in chunk_nums:
			if data is None:
				data = self.chunk(num)
			else:
				data = xor(data, self.chunk(num))
		aDroplet.data = data

		return aDroplet

