#!/usr/bin/python3


from CS312Graph import *
import time


class NetworkRoutingSolver:
    def __init__(self):
        pass

    def initializeNetwork(self, network):
        assert (type(network) == CS312Graph)
        self.network = network

    def getShortestPath(self, destIndex):
        self.dest = destIndex
        path_edges = []
        total_length = 0
        start_node = self.network.nodes[self.source]
        shortest_paths = self.all_paths
        notFound = True
        current = None
        for node in self.network.nodes:
            if node.node_id == destIndex:
                current = node
                break
        if current is None:
            return {'cost': float('inf'), 'path': path_edges}
        previous = None
        while notFound:
            for i in shortest_paths:
                if current.node_id == i.id:
                    total_length += i.dist
                    for node in self.network.nodes:
                        if node.node_id == i.prev:
                            previous = node
                    if previous is None:
                        return {'cost': float('inf'), 'path': path_edges}
                    newEdge = None
                    for edge in previous.neighbors:
                        if edge.src == previous and edge.dest == current:
                            newEdge = edge
                    path_edges.append((current.loc, previous.loc, '{:.0f}'.format(newEdge.length)))
                    if previous == start_node:
                        notFound = False
                    else:
                        current = previous
        return {'cost': total_length, 'path': path_edges}

    # Depending on the implementation we use this function is either O(V) for the array or O(logV) for the heap.
    # Each method is implemented below and has more details for the time complexity.
    # The space complexity is O(V) either way. We store each point and its distance to its connected nodes,
    # along the way we add some variables, like previous, to guide us back to the start, but all space complexity is
    # linear with the amount of points, or V, that we get past in.
    def computeShortestPaths(self, srcIndex, use_heap=False):
        self.source = srcIndex
        t1 = time.time()
        if not use_heap:  # O(V) time complexity (see implementation below)
            distances = []
            unvisited = PriorityArray()
            unvisited.nodes = []
            for index, node in enumerate(self.network.nodes):
                node = new_node(node.node_id, None, float('inf'))
                unvisited.insert(node)
                distances.append(node)
            distances[srcIndex].dist = 0
            unvisited.nodes[srcIndex].dist = 0
            while len(unvisited.nodes) > 0:
                current = unvisited.delete_min()
                graphPoint = self.network.nodes[current.id]
                for i in graphPoint.neighbors:
                    dist = current.dist + i.length
                    if dist < distances[i.dest.node_id].dist:
                        distances[i.dest.node_id].dist = dist
                        distances[i.dest.node_id].prev = i.src.node_id
                        for j in range(len(unvisited.nodes) - 1):
                            if unvisited.nodes[j].id == i.dest.node_id:
                                unvisited.decreaseDist(j, dist)
            self.all_paths = distances
        else:  # Heap implementation: O(logV) time complexity (see implementation below)
            distances = {}
            dist_nodes = []
            for index, node in enumerate(self.network.nodes):
                distances[node.node_id] = float('inf')
                node = new_node(node.node_id, None, float('inf'))
                dist_nodes.append(node)
            distances[srcIndex] = 0
            dist_nodes[srcIndex].dist = 0
            heap = PriorityHeap(distances)
            while len(heap.heap) > 0:
                min_node = heap.delete_min()
                for edge in self.network.nodes[min_node].neighbors:
                    if edge.length + distances[min_node] < distances[edge.dest.node_id]:
                        distances[edge.dest.node_id] = distances[min_node] + edge.length
                        dist_nodes[edge.dest.node_id].dist = distances[min_node] + edge.length
                        dist_nodes[edge.dest.node_id].prev = min_node
                        heap.decreaseKey(edge.dest.node_id, distances[edge.dest.node_id])
            self.all_paths = dist_nodes
        t2 = time.time()
        return t2 - t1


# O(1)
class Node:
    def __init__(self):
        self.id = 0
        self.prev = -1
        self.dist = float('inf')


# O(1)
def new_node(id, prev, dist):
    node = Node()
    node.id = id
    node.prev = prev
    node.dist = dist
    return node


# Dijkstra's Array implementation O(V)
class PriorityArray:
    nodes = []

    # O(1) time complexity since we only call .append
    def insert(self, node):
        print("In A Insert")
        self.nodes.append(node)

    # O(V) time complexity since we have to loop through all our nodes to find the smallest distance
    def delete_min(self):
        print("In A Delete")
        min = float('inf')
        minNode = self.nodes[0]
        for i in self.nodes:
            if i.dist < min:
                min = i.dist
                minNode = i
        self.nodes.remove(minNode)
        return minNode

    # O(1) time complexity since we only change a value that is sent to this function
    def decreaseDist(self, index, newDist):
        print("In A Decrease")
        self.nodes[index].dist = newDist


# Dijkstra's Heap implementation O(logV)
class PriorityHeap:

    def __init__(self, nodes):
        self.heap = []
        self.distances = []
        self.min_index = 0
        self.index_array = list(nodes.keys())
        for node_id, dist in nodes.items():
            self.insert(node_id, dist)

    # O(logV) with the call to bubble up (see bubble up), every other line is constant
    def insert(self, node_id, dist):
        print("In H Insert")
        self.heap.append(node_id)
        self.distances.append(dist)
        self.index_array[node_id] = len(self.heap) - 1
        self.bubble_up(len(self.heap) - 1)

    # O(logV) time complexity since we take the newly added node and only have to look at the nodes
    # directly up to the root node.
    def bubble_up(self, index):
        print("In H Bubble Up")
        sorting = True
        while sorting and index != 0:
            if self.distances[index] < self.distances[(index - 1) // 2]:
                self.switch(index, (index - 1) // 2)
                index = (index - 1) // 2
            else:
                sorting = False

    # O(logV) time complexity since we call bubble_down (see bubble down), every other line is constant
    def delete_min(self):
        print("In H Delete")
        min = self.heap[self.min_index]
        if len(self.heap) == 1:
            del self.heap[0]
            return min
        self.index_array[self.heap[-1]] = 0
        self.index_array[self.heap[0]] = None
        last = self.heap[-1]
        last_dist = self.distances[-1]
        self.heap[self.min_index] = last
        self.distances[self.min_index] = last_dist
        del self.heap[-1]
        del self.distances[-1]
        self.bubble_down(self.min_index)
        return min

    # O(logV) time complexity since we move the last node to the top and only have to look at a single child at each
    # level until there are no more that are smaller, everything inside the loop is constant,
    # including the switch function
    def bubble_down(self, index):
        print("In H Bubble Down")
        while index < (len(self.heap) // 2):  # O(logV) times through
            left_child = float('inf')
            right_child = float('inf')
            if ((2 * index) + 1) < len(self.distances):
                left_child = self.distances[(2 * index) + 1]
            if ((2 * index) + 2) < len(self.distances):
                right_child = self.distances[(2 * index) + 2]
            if left_child <= right_child:
                min_child = left_child
                child_index = (2 * index) + 1
            else:
                min_child = right_child
                child_index = (2 * index) + 2
            if min_child < self.distances[index]:
                self.switch(child_index, index)
                index = child_index
            else:
                return

    # O(logV) since we make a call to bubble up (see bubble up) every other line is constant
    def decreaseKey(self, index, dist):
        print("In H Decrease")
        node = self.index_array[index]
        if node is None:
            return
        self.distances[node] = dist
        self.bubble_up(node)  # O(logV)

    # O(1) since all we do is change the values in our index_array, distances array, and heap
    def switch(self, child_index, parent_index):
        print("In H Switch")
        parent = self.heap[parent_index]
        child = self.heap[child_index]
        self.index_array[parent] = child_index
        self.index_array[child] = parent_index
        temp = self.heap[child_index]
        self.heap[child_index] = self.heap[parent_index]
        self.heap[parent_index] = temp
        temp = self.distances[child_index]
        self.distances[child_index] = self.distances[parent_index]
        self.distances[parent_index] = temp
