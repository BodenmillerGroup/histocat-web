import {v4} from 'uuid';

/**
 * A loader ancestor class containing a default constructor
 * and a stub for the required load() method.
 */
export default class AbstractLoader {
  constructor({
    type, url, requestInit, options,
  }) {
    this.type = type;
    this.url = url;
    this.requestInit = requestInit;
    this.options = options;

    this.subscriptions = {};
  }

  // eslint-disable-next-line class-methods-use-this
  load() {
    throw new Error('The load() method has not been implemented.');
  }

  subscribe(subscriber) {
    const token = v4();
    this.subscriptions[token] = subscriber;
    return token;
  }

  unsubscribe(token) {
    delete this.subscriptions[token];
  }

  publish(data) {
    Object.values(this.subscriptions).forEach((subscriber) => {
      subscriber(data);
    });
  }
}
