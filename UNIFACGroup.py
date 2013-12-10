#!/usr/bin/python


class UNIFACGroup:

    def __init__(self, name, id):
        self.name = name
        self.R = 0.0
        self.Q = 0.0
        self.id = id

    def set_R(self, R):
        self.R = R

    def set_Q(self, Q):
        self.Q = Q

    def get_R(self):
        return self.R

    def get_Q(self):
        return self.Q

    def get_name(self):
        return self.name

    def get_id(self):
        return self.id
