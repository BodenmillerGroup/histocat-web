<template>
  <div class="fill">
    <v-app-bar app dense dark color="primary" clipped-left extension-height="0">
      <v-toolbar-title @click.stop="$router.push({ name: 'groups' }, () => {})" class="toolbar-title">
        {{ appName }}
      </v-toolbar-title>

      <v-spacer />

      <v-btn
        text
        v-if="activeGroupId"
        @click="$router.push({ name: 'group-projects', params: { groupId: activeGroupId } }, () => {})"
      >
        <v-icon left>mdi-alpha-p-box-outline</v-icon>
        Projects
      </v-btn>

      <v-btn
        text
        v-if="activeGroupId && !activeProjectId"
        :to="{ name: 'group-members', params: { groupId: activeGroupId } }"
      >
        <v-icon left>mdi-account-multiple-outline</v-icon>
        Members
      </v-btn>

      <v-btn text v-if="!activeGroupId && !activeProjectId && isAdmin" :to="{ name: 'admin-users' }">
        <v-icon left>mdi-account-outline</v-icon>
        Users
      </v-btn>

      <v-btn text v-if="!activeGroupId && !activeProjectId && isAdmin" :to="{ name: 'admin-models' }">
        <v-icon left>mdi-brain</v-icon>
        Models
      </v-btn>

      <v-spacer />

      <v-menu offset-y tile open-on-hover v-if="activeGroupId && activeProjectId">
        <template v-slot:activator="{ on }">
          <v-btn v-on="on" text small>
            <v-icon left small>mdi-view-dashboard-outline</v-icon>
            Layout
          </v-btn>
        </template>
        <v-list dense>
          <v-list-item v-for="layout in layouts" :key="layout.name" @click="loadLayout(layout.name)">
            <v-list-item-title>{{ layout.name }}</v-list-item-title>
          </v-list-item>
        </v-list>
      </v-menu>

      <v-tooltip bottom v-if="activeGroupId && activeProjectId">
        <template v-slot:activator="{ on }">
          <v-btn v-on="on" icon @click="saveLayout">
            <v-icon>mdi-shape-square-plus</v-icon>
          </v-btn>
        </template>
        <span>Save layout</span>
      </v-tooltip>

      <v-tooltip bottom v-if="activeGroupId && activeProjectId">
        <template v-slot:activator="{ on }">
          <v-btn v-on="on" icon @click="resetLayouts">
            <v-icon>mdi-square-off-outline</v-icon>
          </v-btn>
        </template>
        <span>Reset layouts</span>
      </v-tooltip>

      <v-divider vertical />

      <v-tooltip bottom>
        <template v-slot:activator="{ on }">
          <v-btn v-on="on" icon href="https://plankter.gitbook.io/histocat" target="_blank">
            <v-icon>mdi-help-circle-outline</v-icon>
          </v-btn>
        </template>
        <span>Help</span>
      </v-tooltip>
      <v-menu bottom left offset-y>
        <template v-slot:activator="{ on }">
          <v-btn icon v-on="on">
            <v-icon>mdi-dots-vertical</v-icon>
          </v-btn>
        </template>
        <v-list>
          <v-list-item to="/main/profile">
            <v-list-item-title>Profile</v-list-item-title>
            <v-list-item-action>
              <v-icon>mdi-account</v-icon>
            </v-list-item-action>
          </v-list-item>
          <v-list-item @click="logout">
            <v-list-item-title>Logout</v-list-item-title>
            <v-list-item-action>
              <v-icon>mdi-logout-variant</v-icon>
            </v-list-item-action>
          </v-list-item>
        </v-list>
      </v-menu>
      <ToolbarProgressBar
        :processing="processing"
        :progress="processingProgress"
        :indeterminate="false"
        color="light-blue lighten-2"
        slot="extension"
      />
    </v-app-bar>
    <v-main class="fill">
      <router-view />
    </v-main>
  </div>
</template>

<script lang="ts">
import ToolbarProgressBar from "@/components/ToolbarProgressBar.vue";
import { appName } from "@/env";
import { mainModule } from "@/modules/main";
import { BroadcastManager } from "@/utils/BroadcastManager";
import { WebSocketManager } from "@/utils/WebSocketManager";
import { Component, Vue } from "vue-property-decorator";
import { groupModule } from "@/modules/group";
import { projectsModule } from "@/modules/projects";
import { uiModule } from "@/modules/ui";

const routeGuardMain = async (to, from, next) => {
  if (to.path === "/main") {
    next("/main/groups");
  } else {
    next();
  }
};

@Component({
  components: { ToolbarProgressBar },
})
export default class Main extends Vue {
  readonly uiContext = uiModule.context(this.$store);
  readonly mainContext = mainModule.context(this.$store);
  readonly groupContext = groupModule.context(this.$store);
  readonly projectsContext = projectsModule.context(this.$store);

  readonly appName = appName;

  beforeRouteEnter(to, from, next) {
    routeGuardMain(to, from, next);
  }

  beforeRouteUpdate(to, from, next) {
    routeGuardMain(to, from, next);
  }

  get isAdmin() {
    return this.mainContext.getters.isAdmin;
  }

  get activeGroupId() {
    return this.groupContext.getters.activeGroupId;
  }

  get activeProjectId() {
    return this.projectsContext.getters.activeProjectId;
  }

  async logout() {
    await this.mainContext.actions.userLogOut();
  }

  get processing() {
    return this.uiContext.getters.processing;
  }

  get processingProgress() {
    return this.uiContext.getters.processingProgress;
  }

  get layouts() {
    return this.uiContext.getters.layouts;
  }

  loadLayout(name: string) {
    this.uiContext.mutations.loadLayout(name);
  }

  saveLayout() {
    const name = self.prompt("Please enter layout name:");
    if (name) {
      this.uiContext.mutations.addLayout(name);
    }
  }

  resetLayouts(name: string) {
    this.uiContext.mutations.resetLayouts();
  }

  mounted() {
    WebSocketManager.init(this.$store);
    BroadcastManager.init(this.$store);
  }

  beforeDestroy() {
    WebSocketManager.close();
    BroadcastManager.clear();
  }
}
</script>

<style scoped>
.toolbar-title {
  cursor: pointer;
}
.fill {
  width: 100%;
  height: 100%;
}
</style>

<style>
.link {
  text-decoration: none;
}
</style>
