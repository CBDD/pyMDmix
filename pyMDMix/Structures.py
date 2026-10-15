# To change this template, choose Tools | Templates
# and open the template in the editor.

import os
import time
import errno
import socket

class FileLockException(Exception):
    pass

def _pid_alive(pid):
    "True if a process with this pid exists on this host (or if it cannot be determined)."
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True      # exists but belongs to another user
    except OSError:
        return True
    return True

class FileLock(object):
    """ A file locking mechanism that has context-manager support so
        you can use it in a with statement. This should be relatively cross
        compatible as it doesn't rely on msvcrt or fcntl for the locking.

        The lock file is created next to the locked file and stores the hostname and pid of the owner.
        A lock left behind by a process that died (kill, out-of-memory, crash) is detected and removed:
        on the same host when the owner pid no longer exists, and on any host when the lock is older
        than `stale_after` seconds (no pyMDMix write operation lasts that long).
    """

    def __init__(self, file_name, timeout=10, delay=.05, stale_after=600):
        """ Prepare the file locker. Specify the file to lock and optionally
            the maximum timeout and the delay between each attempt to lock.
        """
        self.is_locked = False
        # Lock next to the target file (absolute path), independent of the current working directory
        self.lockfile = "%s.lock" % os.path.abspath(file_name)
        self.file_name = file_name
        self.timeout = timeout
        self.delay = delay
        self.stale_after = stale_after
        self.fd = None
        self._owner_pid = None

    def _is_stale(self):
        "Check whether the existing lock file belongs to a dead process or is too old."
        try:
            age = time.time() - os.path.getmtime(self.lockfile)
            with open(self.lockfile) as f:
                content = f.read().split()
        except (OSError, ValueError):
            return False         # removed meanwhile or unreadable: let the loop retry
        if len(content) >= 2:
            host, pid = content[0], content[1]
            if host == socket.gethostname() and pid.isdigit() and not _pid_alive(int(pid)):
                return True
        return age > self.stale_after

    def acquire(self):
        """ Acquire the lock, if possible. If the lock is in use, it check again
            every `wait` seconds. It does this until it either gets the lock or
            exceeds `timeout` number of seconds, in which case it throws
            an exception.
        """
        start_time = time.time()
        while True:
            try:
                self.fd = os.open(self.lockfile, os.O_CREAT|os.O_EXCL|os.O_RDWR)
                os.write(self.fd, ("%s %d\n" % (socket.gethostname(), os.getpid())).encode())
                break
            except OSError as e:
                if e.errno != errno.EEXIST:
                    raise
                if self._is_stale():
                    try:
                        os.unlink(self.lockfile)
                    except OSError:
                        pass
                    continue
                if (time.time() - start_time) >= self.timeout:
                    raise FileLockException(
                        "Timeout waiting for lock file %s. Another pyMDMix process is using %s. "
                        "If no other mdmix command is running, the lock is left over from a crashed "
                        "run: delete that .lock file." % (self.lockfile, self.file_name))
                time.sleep(self.delay)
        self.is_locked = True
        self._owner_pid = os.getpid()


    def release(self):
        """ Get rid of the lock by deleting the lockfile.
            When working in a `with` statement, this gets automatically
            called at the end.
        """
        # A forked child (multiprocessing worker) inherits a copy of this object but must never
        # remove the lock of its parent.
        if self.is_locked and self._owner_pid == os.getpid():
            os.close(self.fd)
            try:
                os.unlink(self.lockfile)
            except OSError:
                pass
            self.is_locked = False


    def __enter__(self):
        """ Activated when used in the with statement.
            Should automatically acquire a lock to be used in the with block.
        """
        if not self.is_locked:
            self.acquire()
        return self


    def __exit__(self, type, value, traceback):
        """ Activated at the end of the with statement.
            It automatically releases the lock if it isn't locked.
        """
        if self.is_locked:
            self.release()


    def __del__(self):
        """ Make sure that the FileLock instance doesn't leave a lockfile
            lying around.
        """
        try:
            self.release()
        except Exception:
            pass
