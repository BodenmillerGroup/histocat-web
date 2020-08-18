<template>
  <LoadingView v-if="!group" text="Loading..." />
  <v-container v-else fluid class="px-1 py-0">
    <router-view />
  </v-container>
</template>

<script lang="ts">
import LoadingView from "@/components/LoadingView.vue";
import { WebSocketManager } from "@/utils/WebSocketManager";
import { Component, Vue } from "vue-property-decorator";
import { groupModule } from "@/modules/group";

@Component({
  components: {
    LoadingView,
  },
})
export default class GroupView extends Vue {
  readonly groupContext = groupModule.context(this.$store);

  get group() {
    return this.groupContext.getters.activeGroup;
  }

  async mounted() {
    const groupId = parseInt(this.$router.currentRoute.params.groupId, 10);
    this.groupContext.mutations.setActiveGroupId(groupId);
    await Promise.all([this.groupContext.actions.getMyMember(groupId), this.groupContext.actions.getGroup(groupId)]);
    WebSocketManager.connect(groupId);
  }

  beforeDestroy() {
    WebSocketManager.close();
    this.$store.dispatch("reset");
  }
}
</script>
