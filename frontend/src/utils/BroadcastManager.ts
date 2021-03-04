import { Store } from "vuex";
import { createPubSub, globalPubSub } from "pub-sub-es";

export class BroadcastManager {
  static pubSub: any;

  static init(store: Store<any>) {
    // Creates a new pub-sub event listener stack.
    BroadcastManager.pubSub = createPubSub({ async: false, caseInsensitive: false });
  }

  static publish(event: string, payload?: any, isGlobal: boolean = true) {
    if (isGlobal) {
      globalPubSub.publish(event, payload);
    } else {
      BroadcastManager.pubSub.publish(event, payload);
    }
  }

  static subscribe(event: string, handler: Function, times: number = Infinity, isGlobal: boolean = true) {
    return isGlobal
      ? globalPubSub.subscribe(event, handler, times)
      : BroadcastManager.pubSub.subscribe(event, handler, times);
  }

  static unsubscribe(event: string | object, handler?: Function, isGlobal: boolean = true) {
    if (isGlobal) {
      globalPubSub.unsubscribe(event, handler);
    } else {
      BroadcastManager.pubSub.unsubscribe(event, handler);
    }
  }

  static clear() {
    BroadcastManager.pubSub.clear();
  }
}
