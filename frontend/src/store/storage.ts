import localForage from 'localforage';

const store = localForage.createInstance({
  name: 'histocat',
});

const mapStateForCache = (state) => {
  // Your business logic here
};

export const setState = (state) => {
  const cachedState = mapStateForCache(state);
  return store.setItem('state', cachedState);
};

export const getState = () => store.getItem('state');

export const deleteState = () => store.removeItem('state');
