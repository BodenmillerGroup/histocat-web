<template>
  <router-view></router-view>
</template>

<script lang="ts">
  import { Component, Vue } from 'vue-property-decorator';
  import { store } from '@/store';
  import { readHasAdminAccess } from '@/modules/main/getters';

  const routeGuardAdmin = async (to, from, next) => {
    if (!readHasAdminAccess(store)) {
      next('/main');
    } else {
      next();
    }
  };

  @Component
  export default class Admin extends Vue {
    beforeRouteEnter(to, from, next) {
      routeGuardAdmin(to, from, next);
    }

    beforeRouteUpdate(to, from, next) {
      routeGuardAdmin(to, from, next);
    }
  }
</script>
