<template>
  <v-app id="app">
    <v-content v-if="loggedIn === null">
      <v-container fill-height>
        <v-row align="center" justify="center">
          <v-col>
            <div class="text-center">
              <div class="headline my-12">Loading...</div>
              <v-progress-circular size="100" indeterminate color="primary"></v-progress-circular>
            </div>
          </v-col>
        </v-row>
      </v-container>
    </v-content>
    <router-view v-else />
    <NotificationsManager />
  </v-app>
</template>

<script lang="ts">
import NotificationsManager from "@/components/NotificationsManager.vue";
import { mainModule } from "@/modules/main";
import { Component, Vue } from "vue-property-decorator";

@Component({
  components: {
    NotificationsManager,
  },
})
export default class App extends Vue {
  readonly mainContext = mainModule.context(this.$store);

  get loggedIn() {
    return this.mainContext.getters.isLoggedIn;
  }

  async created() {
    await this.mainContext.actions.checkLoggedIn();
  }
}
</script>
