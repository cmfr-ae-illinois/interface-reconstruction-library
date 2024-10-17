#include "examples/2d_advector/functions.h"

namespace ASHISH {

void P_LinkedList::append(Point value) {
  CoordArray* newNode = new CoordArray(value);
  if (head == nullptr) {
    head = newNode;
    tail = newNode;
    newNode->next = head;  // Point the new node to itself
  } else {
    tail->next = newNode;
    tail = newNode;
    tail->next = head;  // Point the last node to the head
  }
}

void P_LinkedList::setHead(Point value) {
  CoordArray* newNode = new CoordArray(value);
  if (!head) {
    // If the list is empty, set the node as the head and tail, and make it
    // circular.
    head = newNode;
    tail = newNode;
    newNode->next = newNode;
  } else {
    // Insert the node at the beginning of the list.
    newNode->next = head;
    tail->next = newNode;
    head = newNode;
  }
}

void P_LinkedList::setTail(Point value) {
  CoordArray* newNode = new CoordArray(value);
  if (!tail) {
    // If the list is empty, set the node as the head and tail, and make it
    // circular.
    head = newNode;
    tail = newNode;
    newNode->next = newNode;
  } else {
    // Insert the node at the end of the list.
    newNode->next = head;
    tail->next = newNode;
    tail = newNode;
  }
}

CoordArray* P_LinkedList::getHead() { return head; }
CoordArray* P_LinkedList::getTail() { return tail; }

void P_LinkedList::deleteLastNode() {
  // List is empty
  if (!head) {
    return;
  }

  // List has only one node
  if (head == tail) {
    delete head;
    head = nullptr;
    tail = nullptr;
    return;
  }

  // Traverse to find the second last node
  CoordArray* current = head;
  while (current->next != tail) {
    current = current->next;
  }

  // current is now the second last node
  CoordArray* lastNode = tail;
  tail = current;
  tail->next = head;  // Update tail's next to point to head
  delete lastNode;    // Delete the last node
}

void P_LinkedList::printList_t() const {
  if (!head) {
    std::cout << "The list is empty." << std::endl;
    return;
  }
  CoordArray* temp = head;
  do {
    std::cout << "(" << temp->data.x << ", " << temp->data.y << ") " << "-> ";
    temp = temp->next;
  } while (temp != head);
  std::cout << "(head)" << std::endl;
}

void PTS_LinkedList::append(Point_T_sides value) {
  HalfEdge* newNode = new HalfEdge(value);
  if (head == nullptr) {
    head = newNode;
    tail = newNode;
    newNode->next = head;  // Point the new node to itself
  } else {
    tail->next = newNode;
    tail = newNode;
    tail->next = head;  // Point the last node to the head
  }
}

void PTS_LinkedList::setHead(Point_T_sides value) {
  HalfEdge* newNode = new HalfEdge(value);
  if (!head) {
    // If the list is empty, set the node as the head and tail, and make it
    // circular.
    head = newNode;
    tail = newNode;
    newNode->next = newNode;
  } else {
    // Insert the node at the beginning of the list.
    newNode->next = head;
    tail->next = newNode;
    head = newNode;
  }
}

void PTS_LinkedList::setTail(Point_T_sides value) {
  HalfEdge* newNode = new HalfEdge(value);
  if (!tail) {
    // If the list is empty, set the node as the head and tail, and make it
    // circular.
    head = newNode;
    tail = newNode;
    newNode->next = newNode;
  } else {
    // Insert the node at the end of the list.
    newNode->next = head;
    tail->next = newNode;
    tail = newNode;
  }
}

HalfEdge* PTS_LinkedList::getHead() { return head; }
HalfEdge* PTS_LinkedList::getTail() { return tail; }

void PTS_LinkedList::deleteLastNode() {
  if (!head) {
    // List is empty
    return;
  }

  if (head == tail) {
    // List has only one node
    delete head;
    head = nullptr;
    tail = nullptr;
    return;
  }

  // Traverse to find the second last node
  HalfEdge* current = head;
  while (current->next != tail) {
    current = current->next;
  }

  // current is now the second last node
  HalfEdge* lastNode = tail;
  tail = current;
  tail->next = head;  // Update tail's next to point to head
  delete lastNode;    // Delete the last node
}

void PTS_LinkedList::shiftNodesBack() {
  // List is empty or has only one node, no rotation needed
  if (!head || head == tail) {
    return;
  }

  // current =  node just before the current tail
  HalfEdge* current = head;
  while (current->next != tail) {
    current = current->next;
  }

  HalfEdge* newTail = current;
  HalfEdge* newHead = tail;

  // Update pointers
  tail = newTail;
  head = newHead;
  tail->next = head;
}

void PTS_LinkedList::printList_t() const {
  if (!head) {
    std::cout << "The list is empty." << std::endl;
    return;
  }
  HalfEdge* temp = head;
  do {
    std::cout << "(" << temp->data.P.x << ", " << temp->data.P.y << ", "
              << temp->data.t << ", " << temp->data.side_num << ") " << "-> ";
    temp = temp->next;
  } while (temp != head);
  std::cout << "(head)" << std::endl;
}

PTS_LinkedList::~PTS_LinkedList() {
  tail->next = nullptr;
  HalfEdge* temp;
  while (head != nullptr) {
    temp = head;
    head = head->next;
    delete temp;
  }
}

}  // namespace ASHISH