#pragma once

#include <queue>
#include <functional>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <iostream>
#include <stdexcept>

namespace Erwin
{

using namespace std;

template <typename Container>
struct AsyncLoader {

    AsyncLoader( size_t queue_size, function<Container()> f )
        : find_next( f ), size( queue_size ),
          worker_thread( &AsyncLoader<Container>::worker_thread_function, this )
    {
        // first, build the basics
        unique_lock<mutex> lk( m );
        cv_init.wait( lk, [this] { return this->qued.size() >= this->size; } );
    }

    AsyncLoader( const AsyncLoader& ) = delete;
    AsyncLoader( AsyncLoader&& rhs )
        : size( rhs.size ), qued( move( rhs.qued ) ),
          find_next( move( rhs.find_next ) ),
          worker_thread( move( rhs.worker_thread ) ), m(), cv(), cv_init()
    {
    }

    void worker_thread_function()
    {
        while ( true ) {
            unique_lock<mutex> lk( m );
            cv.wait( lk, [this] { return this->qued.size() < this->size; } );
            while ( this->qued.size() < this->size ) {
                try {
                    qued.emplace( find_next() );
                } catch ( domain_error ) {
                    this->size--;
                }
            }
            lk.unlock();
            cv_init.notify_one();
        }
    }


    Container pop_front()
    {
        unique_lock<mutex> lk( m );
        cv_init.wait( lk, [this] { return this->qued.size() >= this->size; } );
        if ( this->qued.empty() == 0 ) throw out_of_range( "out of range" );
        auto front = move( qued.front() );
        qued.pop();
        lk.unlock();
        cv.notify_one();
        return front;
    }

  private:
    size_t size;
    queue<Container> qued;
    function<Container()> find_next;
    thread worker_thread;
    mutex m;
    condition_variable cv;
    condition_variable cv_init;
};

template <typename Container>
AsyncLoader<Container> make_async_loader( size_t queue_size,
                                          function<Container()>&& f )
{
    return AsyncLoader<Container>( queue_size,
                                   forward<function<Container()>>( f ) );
}
}
